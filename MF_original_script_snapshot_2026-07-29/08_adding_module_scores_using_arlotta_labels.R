
# Bootstrap dependencies for a fresh R installation. The package lists are
# explicit so this script is runnable locally or on a cluster without relying
# on a previously prepared library.
cran_packages <- c(
  "dplyr", "Seurat", "viridis", "ggplot2", "stringr", "reshape2",
  "readxl", "rstatix", "data.table", "msigdbr", "ggpubr"
)
bioc_packages <- c("glmGamPoi", "DESeq2")
options(repos = c(CRAN = "https://cloud.r-project.org"))
user_library <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_library)) {
  user_library <- file.path(path.expand("~"), "Documents", "R", "win-library", paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][[1]], sep = "."))
}
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_library, .libPaths())))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_library)
}
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = user_library, dependencies = TRUE)
}
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, lib = user_library, ask = FALSE, update = FALSE)
}
if (!requireNamespace("stats4", quietly = TRUE)) {
  stop("Base R package 'stats4' is unavailable. Repair/reinstall R first.")
}

library(BiocManager)
library(glmGamPoi)
library(dplyr)
library(Seurat)
library(viridis)
library(ggplot2)
library(stringr)
library(reshape2)
library(DESeq2)
library(readxl)
library(rstatix)
library(data.table)
library(msigdbr)
library(ggpubr)


##################################################
###################### The loop ###############
####################################################

celltype_colors <- c(
  "Astroglia" = "#fb9a99",        
  "CFuPN" = "#bebada",           
  "CPN" = "#a50f15",             
  "Glial precursors" = "#8dd3c7",
  "IN progenitors" = "#fccde5",  
  "IP" = "#02818a",             
  "Immature IN" = "#dfc27d",     
  "Unspecified PN" = "#8c510a",  
  "aRG" = "#08519C",             
  "oRG" = "#f6e8c3",             
  "oRG/Astroglia" = "#c7e9c0"   
)

# Testing with static included also
setwd("~/Microfluidics McCain Collab/single_cell_run")
seur<- readRDS('Final_joined_obj_label_transfer_Arlotta.Robj')

seur
head(seur)
tail(seur)

# Load Seurat object
setwd("J:/MF_Single_Cell_Updated_No_Static")
seur <- readRDS("Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj")
expressed_genes <- rownames(seur)

unique(seur@meta.data$predicted_CellType)

#  Load Hallmark gene sets
target_gs_names <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS"
)

# Now filter 
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(gs_name %in% target_gs_names)


#  Filter gene sets to only include expressed genes
hallmark_filtered <- hallmark_df %>%
  filter(gene_symbol %in% expressed_genes)

#  Split into list and drop empty sets
hallmark_list <- split(hallmark_filtered$gene_symbol, hallmark_filtered$gs_name)
hallmark_list <- hallmark_list[lengths(hallmark_list) > 0]  # Drop empty

# Clean names for scoring
names(hallmark_list) <- gsub("HALLMARK_", "", names(hallmark_list))

# Add module scores
seur <- AddModuleScore(seur, features = hallmark_list, name = "Score")

# Fix duplicate column names
score_cols <- grep("Score[0-9]+$", colnames(seur@meta.data), value = TRUE)
names(score_cols) <- names(hallmark_list)
colnames(seur@meta.data)[colnames(seur@meta.data) %in% score_cols] <- paste0("Score_", names(score_cols))

# Use updated score names
module_scores <- paste0("Score_", names(score_cols))

# Metadata and plotting setup
meta <- seur@meta.data
cell_types <- c("CPN", "CFuPN", 'oRg/Astroglia', 'Unspecified PN', 'aRG', 'IP', 'oRG', 'Astroglia', 'Glial precursors')
setwd("J:/MF_Single_Cell_Updated_No_Static")
output_dir <- "Module_Score_Plots_LM_By_CellType_wStatic_wArlotta_5mo"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Main loop
for (score in module_scores) {
  for (ct in cell_types) {
    
    # Step 1: Average score per sample × cell type
    module_avg <- meta %>%
      group_by(sample, predicted_CellType) %>%
      summarise(Average_Score = mean(!!sym(score), na.rm = TRUE), .groups = "drop") %>%
      filter(predicted_CellType == ct)
    
    # Step 2: Assign Batch_Group
    df <- module_avg %>%
      mutate(Batch_Group = case_when(
        grepl("microfluidic", sample, ignore.case = TRUE) ~ "Microfluidic",
        grepl("bioreactor", sample, ignore.case = TRUE) ~ "Bioreactor"
      )) %>%
      filter(!is.na(Batch_Group)) %>%
      mutate(Batch_Group = factor(Batch_Group, levels = c("Bioreactor", "Microfluidic")))
    
    # Step 3: Skip if insufficient data
    if (nrow(df) < 4 || length(unique(df$Batch_Group)) < 2) {
      message(paste("Skipping", score, ct, "- insufficient data"))
      next
    }
    
    # Step 4: Run linear model
    model <- lm(Average_Score ~ Batch_Group, data = df)
    model_summary <- summary(model)
    pval <- tryCatch(coef(model_summary)["Batch_GroupMicrofluidic", "Pr(>|t|)"], error = function(e) NA)
    
    # Step 5: Format p-value label
    pval_label <- if (!is.na(pval)) {
      if (pval < 0.001) "***"
      else if (pval < 0.01) "**"
      else if (pval < 0.05) "*"
      else "ns"
    } else {
      "NA"
    }
    
    # Step 6: Create p-value data.frame
    p_df <- data.frame(
      group1 = "Bioreactor",
      group2 = "Microfluidic",
      p = pval,
      y.position = max(df$Average_Score, na.rm = TRUE) * 1.1,
      label = paste0("p = ", signif(pval, 3), " (", pval_label, ")")
    )
    
    
    print(levels(df$Batch_Group))
    print(p_df)
    
    # Step 7: Build plot
    p <- ggplot(df, aes(x = Batch_Group, y = Average_Score, fill = Batch_Group)) +
      geom_boxplot(width = 0.5, color = "black", size = 1, alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 4, alpha = 0.8, color = "black") +
      scale_fill_manual(values = c("Microfluidic" = "#B8E2DE", "Bioreactor" = "#C9DCF4")) +
      theme_classic() +
      theme(
        axis.line = element_line(linewidth = 1, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      labs(
        title = paste(score, "in", ct),
        x = "Condition",
        y = paste("Average", score, "Score")
      ) +
      stat_pvalue_manual(p_df, label = "label", size = 5, inherit.aes = FALSE)
    
    
    
    # Step 8: Save plot
    ggsave(
      filename = file.path(output_dir, paste0(score, "_", ct, "_lm_score_plot.pdf")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}




#####################################################################################################################################

# Load all GO Biological Process terms
#go_bp <- msigdbr(
#  species = "Homo sapiens",
#  collection = "C5", 
#  subcollection = "GO:BP"
#)
# Filter for neuronal maturation-related pathways
#neuronal_terms <- go_bp %>%
#  filter(grepl("maturation|development|differentiation.*neuron|neuron.*differentiation", 
#               gs_name, ignore.case = TRUE))
# View distinct pathway names
#unique(neuronal_terms$gs_name)

#neuronal_gene_sets <- neuronal_terms %>%
#  split(x = .$gene_symbol, f = .$gs_name)

######### Go BP subsetting as above results in a large list of Go terms, too many to look at/ add modules with I think 05 28 25

#####################################################################################################################################

#loop for running comparison on all conditions: bioreactor, static, microfluidic, arlotta



library(dplyr)
library(ggplot2)
library(ggpubr)


# Colors (you gave MF + BR; I’m adding reasonable defaults for Static + Arlotta)
group_colors <- c(
  "Microfluidic" = "#B8E2DE",
  "Bioreactor"   = "#C9DCF4",
  "Static"       = "#1b9e77",
  "Arlotta"      = "#7570b3"
)

group_levels <- c("Static", "Bioreactor", "Microfluidic", "Arlotta")

for (score in module_scores) {
  for (ct in cell_types) {
    
    # Step 1: Average score per sample × cell type
    module_avg <- meta %>%
      group_by(sample, predicted_CellType) %>%
      summarise(Average_Score = mean(.data[[score]], na.rm = TRUE), .groups = "drop") %>%
      filter(predicted_CellType == ct)
    
    # Step 2: Assign Batch_Group (4 groups)
    df <- module_avg %>%
      mutate(
        Batch_Group = case_when(
          grepl("static",      sample, ignore.case = TRUE) ~ "Static",
          grepl("bioreactor",  sample, ignore.case = TRUE) ~ "Bioreactor",
          grepl("microfluidic",sample, ignore.case = TRUE) ~ "Microfluidic",
          grepl("arlotta",     sample, ignore.case = TRUE) ~ "Arlotta",
          TRUE ~ NA_character_
        ),
        Batch_Group = factor(Batch_Group, levels = group_levels)
      ) %>%
      filter(!is.na(Batch_Group))
    
    # Step 3: Skip if insufficient data
    if (nrow(df) < 4 || n_distinct(df$Batch_Group) < 2) {
      message(paste("Skipping", score, ct, "- insufficient data"))
      next
    }
    
    # Step 4: Stats
    # If only 2 groups present -> t-test. If >=3 groups -> ANOVA + Tukey.
    ngroups <- n_distinct(df$Batch_Group)
    
    if (ngroups == 2) {
      present <- levels(droplevels(df$Batch_Group))
      comps <- list(present)
      stat_method <- "t.test"
      
      p_df <- compare_means(
        formula = Average_Score ~ Batch_Group,
        data = df,
        method = stat_method
      ) %>%
        mutate(
          group1 = present[1],
          group2 = present[2],
          y.position = max(df$Average_Score, na.rm = TRUE) * 1.10,
          label = paste0("p = ", signif(p, 3),
                         " (", ifelse(p < 0.001, "***",
                                      ifelse(p < 0.01, "**",
                                             ifelse(p < 0.05, "*", "ns"))), ")")
        ) %>%
        select(group1, group2, y.position, label)
      
    } else {
      # ANOVA
      aov_fit <- aov(Average_Score ~ Batch_Group, data = df)
      # Tukey (pairwise) and format for stat_pvalue_manual
      tuk <- TukeyHSD(aov_fit)$Batch_Group
      tuk_df <- as.data.frame(tuk)
      tuk_df$comparison <- rownames(tuk_df)
      
      # Split "A-B" into group1/group2
      tuk_df <- tuk_df %>%
        tidyr::separate(comparison, into = c("group2", "group1"), sep = "-", remove = FALSE) %>%
        mutate(
          p = `p adj`,
          # stack brackets above the data
          y.position = max(df$Average_Score, na.rm = TRUE) * 1.05 +
            (row_number() * 0.03 * max(df$Average_Score, na.rm = TRUE)),
          label = paste0("p = ", signif(p, 3),
                         " (", ifelse(p < 0.001, "***",
                                      ifelse(p < 0.01, "**",
                                             ifelse(p < 0.05, "*", "ns"))), ")")
        )
      
      # Keep only comparisons where both groups exist in df (robust to missing groups)
      present_groups <- unique(as.character(df$Batch_Group))
      p_df <- tuk_df %>%
        filter(group1 %in% present_groups & group2 %in% present_groups) %>%
        select(group1, group2, y.position, label)
    }
    
    # Step 5: Plot
    p <- ggplot(df, aes(x = Batch_Group, y = Average_Score, fill = Batch_Group)) +
      geom_boxplot(width = 0.55, color = "black", size = 1, alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.18, size = 3.8, alpha = 0.8, color = "black") +
      scale_fill_manual(values = group_colors, drop = TRUE) +
      theme_classic() +
      theme(
        axis.line = element_line(linewidth = 1, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      labs(
        title = paste(score, "in", ct),
        x = "Condition",
        y = paste("Average", score, "Score")
      ) +
      stat_pvalue_manual(
        p_df,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        size = 4,
        inherit.aes = FALSE
      )
    
    # Step 6: Save
    ggsave(
      filename = file.path(output_dir, paste0(score, "_", ct, "_score_plot_4groups.pdf")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}





