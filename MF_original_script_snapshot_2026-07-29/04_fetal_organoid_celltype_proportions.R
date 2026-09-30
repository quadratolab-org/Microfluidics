

########## Cell Type Proportions of Fetal and Organoid
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)




#load('tb.percentage.final.Robj')
#head(tb.final.all)
#head(reference)
#tb.filtered <- tb.final.all %>%
#filter(!grepl("bio|micro|stat", donor.id, ignore.case = TRUE))

#unique(tb.filtered$donor.id)

#setwd("C:/Users/jeanp/Downloads")
#load('object.reference.wang.down.3000.Robj')


setwd("N:/JoJo/Team Microfluidic-JP/seurat.object/Wang.trimester.static.micro.bio.arlotta")


load("object.merge.arlotta.wang.integration.rpca.cluster.Robj")

head(object.rpca)
unique(object.rpca$data.set)
unique(object.rpca$subclass)

object.rpca.sub <- subset(object.rpca, subset = subclass %in% c('Glutamatergic neuron', 'Radial glia', 'IPC-EN', 'IPC-Glia', 'Astrocyte'))

head(object.rpca.sub)

meta_subset <- as.data.frame(object.rpca.sub@meta.data) %>%
  dplyr::select(data.set, Group, subclass, donor_id, type_updated)

head(meta_subset)


meta_pct_donor <- meta_subset %>%
  group_by(Group, donor_id, subclass) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Group, donor_id) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()


head(meta_pct_donor)



setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")
seur2 <- readRDS("Final_joined_obj_label_transfer_Arlotta.Robj")
seur2
#seur2
#An object of class Seurat 
#44206 features across 70846 samples within 1 assay 
#Active assay: RNA (44206 features, 2000 variable features)
#3 layers present: data, counts, scale.data
#3 dimensional reductions calculated: pca, integrated.cca, umap.cca



seur3 <- subset(seur2, subset = data.set == "quadrato" )
seur3
#seur3
#An object of class Seurat 
#44206 features across 39059 samples within 1 assay 
#Active assay: RNA (44206 features, 2000 variable features)
#3 layers present: data, counts, scale.data
#3 dimensional reductions calculated: pca, integrated.cca, umap.cca

head(seur3)
seur4 <- subset(seur3, subset = sample.simple != "Static")
seur4

Idents(seur4)<- 'predicted_CellType'
unique(Idents(seur4))
seur4 <- subset(seur4, subset = predicted_CellType %in% c("CFuPN", "aRG", "CPN", "oRG", "IP", "Glial precursors", "Astroglia", "oRG/Astroglia"))
seur4
head(seur4)



new_groups <- c(
  'CFuPN' = "Glutamatergic neuron",
  'CPN' = "Glutamatergic neuron",
  'aRG' = "Radial Glia",
  'oRG' = "Radial Glia",
  'IP' = "IPC-EN",
  'Glial precursors' = "IPC-Glia",
  'Astroglia' = "Astrocyte",
  'oRG/Astroglia' = "Radial Glia"
)

# Map current identities to new groups
seur4$celltype <- plyr::mapvalues(
  x = as.character(Idents(seur4)),
  from = names(new_groups),
  to = unname(new_groups)
)



Idents(seur4)<- 'celltype'
unique(Idents(seur4))

seur4
head(seur4)




cell_names <- Cells(seur4)  

# 2. Downsample per sample group
min_cells <- min(table(seur4$sample))
set.seed(123)

cells_to_keep <- seur4@meta.data %>%
  mutate(cell_id = Cells(seur4)) %>%  
  group_by(sample) %>%
  sample_n(min_cells) %>%
  pull(cell_id)


seur4_downsampled <- subset(seur4, cells = cells_to_keep)


# Compute percentages per donor.id and celltype
meta_df <- as.data.frame(seur4_downsampled@meta.data)

seur4_pct <- meta_df %>%
  dplyr::count(sample, celltype) %>%
  group_by(sample) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

head(seur4_pct)
head(meta_pct_donor)

meta_pct_donor <- meta_pct_donor %>%
  dplyr::rename(
    celltype = subclass,
    sample = donor_id
  )


# Label and clean filtered table
meta_pct_donor_labeled <- meta_pct_donor %>%
  dplyr::mutate(source = "fetal") %>%
  dplyr::select(Group, sample, celltype, n, percentage, source)

head(meta_pct_donor_labeled)

head(seur4_pct)
# Label and clean Seurat percentages
seur_labeled <- seur4_pct %>%
  dplyr::mutate(source = "organoid") %>%
  dplyr::select(celltype, sample, percentage, source)

head(seur_labeled)

# Combine both
combined_long <- dplyr::bind_rows(meta_pct_donor_labeled, seur_labeled)

head(combined_long)
table(combined_long$source)
class(combined_long)

combined_df <- combined_long %>%
  mutate(celltype = ifelse(celltype == "Radial Glia", "Radial glia", celltype))



#setwd("~/Microfluidics McCain Collab/single_cell_run/Wang_fetal_proportions")
#saveRDS(combined_long, file = "combined_df_fetal_organoid_celltypes_Wang.rds")
#write.csv(combined_long, file = 'combined_df_fetal_organoid_celltypes_Wang.csv' )




# Load your CSV
combined_df <- read_csv("combined_df_fetal_organoid_celltypes_Wang.csv")

# 1. Add organoid_source column using sample string
combined_df <- combined_df %>%
  mutate(
    organoid_source = case_when(
      grepl("^bioreactor", sample, ignore.case = TRUE) ~ "bioreactor",
      grepl("^microfluidic", sample, ignore.case = TRUE) ~ "microfluidic",
      TRUE ~ NA_character_
    )
  )

# 2. Fix cell type naming
combined_df <- combined_df %>%
  mutate(celltype = ifelse(celltype == "Radial Glia", "Radial glia", celltype))

# 3. Filter only fetal stages of interest
combined_df <- combined_df %>%
  filter(!Group %in% c("Infancy", "Adolescence"))

# 4. Set factor levels
combined_df$Group <- factor(combined_df$Group, levels = c("First_trimester", "Second_trimester", "Third_trimester"))
combined_df$organoid_source <- factor(combined_df$organoid_source, levels = c("bioreactor", "microfluidic"))

# Loop over each cell type and generate plots
celltypes <- unique(combined_df$celltype)



combined_df <- combined_df %>%
  mutate(
    organoid_source = case_when(
      grepl("^bioreactor", sample, ignore.case = TRUE) ~ "bioreactor",
      grepl("^microfluidic", sample, ignore.case = TRUE) ~ "microfluidic",
      TRUE ~ NA_character_
    ),
    organoid_source = factor(organoid_source, levels = c("bioreactor", "microfluidic"))
  )


head(combined_df)




for (ct in celltypes) {
  # Subset for current celltype
  plot_data <- combined_df %>% filter(celltype == ct)
  
  # Re-extract organoid data and explicitly assign factors from sample
  organoid_data <- plot_data %>%
    filter(!is.na(organoid_source))
  
   
  
  # Fetal summary: mean + SD per trimester
  fetal_summary <- plot_data %>%
    filter(source == "fetal") %>%
    group_by(Group) %>%
    summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      sd_pct = sd(percentage, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(organoid_data, aes(x = organoid_source, y = percentage, fill = organoid_source)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6, color= NA) +  
    geom_hline(data = fetal_summary, aes(yintercept = mean_pct),
               linetype = "dotted", color = "#0072B2", linewidth = 1) +
    geom_rect(data = fetal_summary,
              aes(xmin = 0.5, xmax = 2.5,
                  ymin = mean_pct - sd_pct,
                  ymax = mean_pct + sd_pct),
              inherit.aes = FALSE,
              fill = "#0072B2", alpha = 0.08) +
    geom_text(data = fetal_summary,
              aes(x = 2.75, y = mean_pct, label = Group),  # shifted further right
              inherit.aes = FALSE,
              hjust = 0,
              size = 4,
              fontface = "italic",
              color = "#0072B2") +
    scale_fill_manual(values = c("bioreactor" = "#E69F00", "microfluidic" = "#9E79B9")) +
    scale_y_continuous(limits = c(0, NA)) +
    coord_cartesian(clip = "off") +  # allow labels to go outside plotting area
    labs(
      title = paste(ct),
      x = NULL,
      y = "Percent"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin = margin(10, 100, 10, 10)  # extended right margin for labels
    )
  
  # Save the plot
  ggsave(
    filename = paste0("celltype_boxplot_", gsub(" ", "_", ct), ".pdf"),
    plot     = p,
    width    = 5,
    height   = 3,
    device   = cairo_pdf   
  )
  
  
}







##############################################################################################################

############ Updated analysis using JoJo's label transfer ####################################################

##############################################################################################################



setwd("N:/JoJo/Team Microfluidic-JP/seurat.object/Quadrato.Arlotta.object.with.prediction")


load("quadrato.only.Robj")


head(quadrato)
tail(quadrato)

quadrato


setwd("N:/JoJo/Team Microfluidic-JP/seurat.object/Wang.trimester.static.micro.bio.arlotta")

load("object.merge.arlotta.wang.integration.rpca.cluster.Robj")

head(object.rpca)
tail(object.rpca)



####################################################################################################################

# ============================================================
# Organoid vs Fetal Trimester Comparison (Replicate-level + cleanup)
# ============================================================

library(Seurat)
library(dplyr)
library(lme4)
library(ggplot2)
library(ggpubr)
library(broom.mixed)
library(purrr)

# Helper -----------------------------------------------------
clean_names <- function(x) {
  x %>%
    as.character() %>%
    gsub(" |/|-", "_", .) %>%
    gsub("[^0-9A-Za-z_]", "", .) %>%
    factor()
}


# Subset fetal-only cells
fetal_only <- subset(object.rpca, subset = !(data.set %in% c("quadrato", "arlotta")))
head(fetal_only)
rm(object.rpca)
# Unified metadata -------------------------------------------
meta_org <- quadrato@meta.data %>%
  mutate(
    replicate = sample,
    ConditionGroup = case_when(
      grepl("Bioreactor", sample, ignore.case = TRUE)   ~ "Bioreactor",
      grepl("Microfluidic", sample, ignore.case = TRUE) ~ "Microfluidic",
      grepl("Static", sample, ignore.case = TRUE)       ~ "Static",
      TRUE ~ NA_character_
    ),
    subclass = predicted.id.subclass
  ) %>%
  filter(!is.na(ConditionGroup)) %>%
  mutate(ConditionGroup = factor(ConditionGroup,
                                 levels = c("Bioreactor","Microfluidic","Static")))

meta_fetal <- fetal_only@meta.data %>%
  mutate(
    replicate = donor_id,
    ConditionGroup = Group %||% Group
  ) %>%
  filter(ConditionGroup %in% c("First_trimester","Second_trimester","Third_trimester")) %>%
  mutate(ConditionGroup = factor(ConditionGroup,
                                 levels = c("First_trimester","Second_trimester","Third_trimester")))

metadata <- bind_rows(meta_org, meta_fetal) %>%
  mutate(subclass = clean_names(subclass))

# Remove unwanted fetal cell types ---------------------------
remove_celltypes <- c("Microglia", "Oligodendrocyte", "OPC", "Unknown", "Vascular")
metadata <- metadata %>%
  filter(!subclass %in% clean_names(remove_celltypes))

# Output directory --------------------------------------------
plots_dir <- "J:/MF and Fetal Data Analysis"
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Pairwise comparisons ---------------------------------------
pairwise_comparisons <- expand.grid(
  Organoid = c("Bioreactor", "Microfluidic", "Static"),
  Trimester = c("First_trimester", "Second_trimester", "Third_trimester"),
  stringsAsFactors = FALSE
)

# GLMM + replicate plotting -----------------------------------
glmm_results_all <- list()

for (i in 1:nrow(pairwise_comparisons)) {
  organoid_group <- pairwise_comparisons$Organoid[i]
  fetal_group    <- pairwise_comparisons$Trimester[i]
  comparison_name <- paste0(organoid_group, "_vs_", fetal_group)
  cat("Running comparison:", comparison_name, "\n")
  
  meta_sub <- metadata %>%
    filter(ConditionGroup %in% c(organoid_group, fetal_group)) %>%
    mutate(ConditionGroup = factor(ConditionGroup, levels = c(fetal_group, organoid_group)))
  
  valid_celltypes <- levels(droplevels(meta_sub$subclass))
  
  # GLMMs -----------------------------------------------------
  glmm_results <- map_dfr(valid_celltypes, function(ct) {
    df_ct <- meta_sub %>%
      mutate(is_ct = as.integer(subclass == ct))
    
    if (n_distinct(df_ct$ConditionGroup) < 2 || n_distinct(df_ct$replicate) < 2) {
      return(tibble(CellType = ct, p_value = NA_real_))
    }
    
    ctl <- glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6),
      check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
      check.conv.grad     = .makeCC(action = "ignore", tol = 1e-3)
    )
    
    full_mod <- tryCatch(
      glmer(is_ct ~ ConditionGroup + (1 | replicate),
            data = df_ct, family = binomial, control = ctl),
      error = function(e) NULL
    )
    null_mod <- tryCatch(
      glmer(is_ct ~ 1 + (1 | replicate),
            data = df_ct, family = binomial, control = ctl),
      error = function(e) NULL
    )
    
    pval <- if (!is.null(full_mod) && !is.null(null_mod)) {
      anova(null_mod, full_mod, test = "Chisq")$`Pr(>Chisq)`[2]
    } else NA_real_
    
    tibble(CellType = ct, p_value = pval)
  }) %>%
    mutate(
      adj_p_value = p.adjust(p_value, method = "BH"),
      Comparison  = comparison_name
    )
  
  glmm_results_all[[comparison_name]] <- glmm_results
  
  # Compute replicate-level proportions ------------------------
  df_prop <- meta_sub %>%
    group_by(replicate, ConditionGroup, subclass) %>%
    summarise(cell_count = n(), .groups = "drop_last") %>%
    mutate(total_cells = sum(cell_count)) %>%
    ungroup() %>%
    mutate(proportion = cell_count / total_cells)
  
  # Define color palette (distinct, publication quality)
  condition_colors <- c(
    "Bioreactor"      = "#D81B60",  # magenta-red
    "Microfluidic"    = "#1E88E5",  # blue
    "Static"          = "#004D40",  # teal-green
    "First_trimester" = "#FFC107",  # amber
    "Second_trimester"= "#F57C00",  # orange
    "Third_trimester" = "#7B1FA2"   # violet
  )
  
  # Plot each replicate ---------------------------------------
  for (ct in unique(df_prop$subclass)) {
    ct_data <- df_prop %>% filter(subclass == ct)
    if (nrow(ct_data) == 0) next
    
    # Group summary (mean ± SD)
    df_summary <- ct_data %>%
      group_by(ConditionGroup) %>%
      summarise(mean_prop = mean(proportion, na.rm = TRUE),
                sd_prop   = sd(proportion, na.rm = TRUE),
                .groups   = "drop")
    
    pval_info <- glmm_results %>% filter(CellType == ct) %>% pull(adj_p_value)
    if (length(pval_info) != 1) pval_info <- NA_real_
    pval_label <- paste0("Adj. p = ", signif(pval_info, 3))
    
    y_max <- max(ct_data$proportion, na.rm = TRUE) * 1.35
    x_center <- mean(seq_along(unique(ct_data$replicate)))
    
    p <- ggplot(ct_data, aes(x = replicate, y = proportion, fill = ConditionGroup)) +
      # individual replicate bars
      geom_col(color = "black", width = 0.75, alpha = 0.9) +
      # overlay group mean ± SD
      geom_errorbar(
        data = df_summary,
        aes(x = ConditionGroup,
            ymin = mean_prop - sd_prop,
            ymax = mean_prop + sd_prop),
        inherit.aes = FALSE,
        width = 0.25, linewidth = 0.9, color = "black"
      ) +
      geom_point(
        data = df_summary,
        aes(x = ConditionGroup, y = mean_prop),
        inherit.aes = FALSE,
        shape = 21, size = 4, fill = "white", color = "black"
      ) +
      scale_fill_manual(values = condition_colors, drop = FALSE) +
      labs(
        title = paste0(ct, " | ", comparison_name),
        x = "Replicate",
        y = "Proportion per replicate"
      ) +
      annotate("text",
               x = x_center,
               y = y_max * 0.92,
               label = pval_label,
               fontface = "bold",
               size = 5) +
      theme_classic(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title  = element_text(size = 12, face = "bold"),
        plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        plot.margin = margin(10, 15, 10, 15)
      ) +
      ylim(0, y_max)
    
    out_file <- file.path(plots_dir,
                          paste0("Prop_", ct, "_", comparison_name, "_replicates.pdf"))
    ggsave(out_file, plot = p, width = 6.8, height = 4.8)
  }
  
}

# Combine & Save Results ------------------------------------
glmm_results_combined <- bind_rows(glmm_results_all)
write.csv(glmm_results_combined,
          file.path(plots_dir, "GLMM_Organoid_vs_Trimester_results.csv"),
          row.names = FALSE)

cat("✅ Analysis complete! Results and plots saved in:", plots_dir, "\n")




###################################################################
##################################################################################
# ============================================
# FACET: Organoid vs All Fetal Trimesters
# ============================================

library(ggsignif)  

facet_dir <- file.path(plots_dir, "Facet_by_Trimester_v2")
if (!dir.exists(facet_dir)) dir.create(facet_dir, recursive = TRUE)

# same color palette as before
condition_colors <- c(
  "Bioreactor"      = "#D81B60",
  "Microfluidic"    = "#1E88E5",
  "Static"          = "#004D40",
  "First_trimester" = "#FFC107",
  "Second_trimester"= "#F57C00",
  "Third_trimester" = "#7B1FA2"
)

for (org in c("Bioreactor", "Microfluidic", "Static")) {
  
  df_all <- metadata %>%
    filter(ConditionGroup %in% c("First_trimester","Second_trimester","Third_trimester", org)) %>%
    mutate(GroupType = ifelse(ConditionGroup %in% c("First_trimester","Second_trimester","Third_trimester"),
                              "Fetal", "Organoid"))
  
  for (ct in unique(df_all$subclass)) {
    df_ct <- df_all %>% filter(subclass == ct)
    
    # Compute proper per-replicate proportions relative to total cells
    df_prop_all <- df_all %>%
      group_by(replicate, ConditionGroup) %>%
      summarise(total_cells = n(), .groups = "drop")
    
    df_prop_ct <- df_ct %>%
      group_by(replicate, ConditionGroup) %>%
      summarise(cell_count = n(), .groups = "drop") %>%
      left_join(df_prop_all, by = c("replicate", "ConditionGroup")) %>%
      mutate(proportion = cell_count / total_cells)
    
    # mean ± sd per group
    df_summary <- df_prop_ct %>%
      group_by(ConditionGroup) %>%
      summarise(mean_prop = mean(proportion, na.rm = TRUE),
                sd_prop   = sd(proportion, na.rm = TRUE),
                .groups = "drop")
    
    # derive p-values for all fetal vs this organoid
    comparison_rows <- glmm_results_combined %>%
      filter(grepl(org, Comparison) & CellType == ct)
    comp_df <- tibble(
      group1 = c("First_trimester","Second_trimester","Third_trimester"),
      group2 = org,
      pval   = comparison_rows$adj_p_value[match(group1, gsub(".*_vs_", "", comparison_rows$Comparison))]
    ) %>%
      filter(!is.na(pval)) %>%
      mutate(y_position = max(df_prop_ct$proportion, na.rm = TRUE) * (1.05 + 0.05 * (row_number() - 1)))
    
    # build plot
    y_max <- max(df_prop_ct$proportion, na.rm = TRUE) * 1.35
    
    p <- ggplot() +
      # transparent bars for each condition
      geom_col(data = df_summary,
               aes(x = ConditionGroup, y = mean_prop, fill = ConditionGroup),
               color = "black", alpha = 0.7, width = 0.7) +
      # replicate points
      geom_jitter(data = df_prop_ct,
                  aes(x = ConditionGroup, y = proportion),
                  color = "black", width = 0.15, size = 2.2, alpha = 0.9) +
      # mean ± sd overlay
      geom_errorbar(data = df_summary,
                    aes(x = ConditionGroup,
                        ymin = mean_prop - sd_prop,
                        ymax = mean_prop + sd_prop),
                    width = 0.25, color = "black", linewidth = 0.8) +
      scale_fill_manual(values = condition_colors, drop = FALSE) +
      ylim(0, y_max) +
      labs(
        title = paste0(ct, " | ", org, " vs. Fetal Trimesters"),
        x = NULL,
        y = "Proportion per replicate"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "none"
      )
    
    # add p-value comparisons (black lines + labels)
    if (nrow(comp_df) > 0) {
      p <- p + geom_signif(
        data = comp_df,
        aes(xmin = group1, xmax = group2, annotations = paste0("Adj. p=", signif(pval, 2)),
            y_position = y_position),
        manual = TRUE, tip_length = 0.01, textsize = 3.8, vjust = -0.3,
        color = "black", size = 0.7
      )
    }
    
    out_file <- file.path(facet_dir, paste0("Facet_", ct, "_", org, "_vs_Trimesters_bars.pdf"))
    ggsave(out_file, plot = p, width = 7.2, height = 5)
  }
}



##########################################################################################################
#############################################################################

# ================================================================
# Hierarchically clustered heatmap of prediction max scores
# ================================================================

library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# --- Extract and prepare metadata --------------------------------

meta <- quadrato@meta.data %>%
  dplyr::select(sample, predicted.id.subclass, prediction.score.max.subclass) %>%
  rename(CellType = predicted.id.subclass,
         PredScore = prediction.score.max.subclass) %>%
  filter(!is.na(CellType), !is.na(PredScore))

# --- Calculate mean max prediction score per cell type × sample ---
heat_df <- meta %>%
  group_by(CellType, sample) %>%
  summarise(mean_score = mean(PredScore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = mean_score, values_fill = 0) %>%
  column_to_rownames("CellType")

# --- Define color palette ----------------------------------------
# e.g. white to blue gradient
heat_colors <- colorRampPalette(brewer.pal(9, "Blues"))(100)

# --- Optional: annotation for sample types -----------------------
sample_ann <- data.frame(
  Condition = case_when(
    grepl("Bioreactor", colnames(heat_df), ignore.case = TRUE)   ~ "Bioreactor",
    grepl("Microfluidic", colnames(heat_df), ignore.case = TRUE) ~ "Microfluidic",
    grepl("Static", colnames(heat_df), ignore.case = TRUE)       ~ "Static",
    TRUE                                                         ~ "Unknown"
  )
)
rownames(sample_ann) <- colnames(heat_df)

ann_colors <- list(
  Condition = c(
    Bioreactor   = "#D81B60",
    Microfluidic = "#1E88E5",
    Static       = "#004D40",
    Unknown      = "grey80"
  )
)

# --- Plot heatmap -------------------------------------------------
pheatmap(
  mat               = as.matrix(heat_df),
  color             = heat_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  scale             = "none",
  border_color      = NA,
  cellwidth         = 18,
  cellheight        = 14,
  annotation_col    = sample_ann,
  annotation_colors = ann_colors,
  fontsize_row      = 10,
  fontsize_col      = 10,
  main              = "Mean Prediction Confidence per Cell Type (Quadrato)"
)

# --- Optional: Save to file --------------------------------------
pdf("J:/MF and Fetal Data Analysis/PredictionScore_Heatmap.pdf",
    width = 9, height = 6)
pheatmap(
  mat               = as.matrix(heat_df),
  color             = heat_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  scale             = "none",
  border_color      = NA,
  cellwidth         = 18,
  cellheight        = 14,
  annotation_col    = sample_ann,
  annotation_colors = ann_colors,
  fontsize_row      = 10,
  fontsize_col      = 10,
  main              = "Mean Prediction Confidence per Cell Type (Quadrato)"
)
dev.off()








