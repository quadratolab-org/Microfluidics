
library(dplyr)
library(lme4)
library(ggplot2)
library(ggpubr)
library(broom.mixed)
library(purrr)

# Helper to clean factor levels to R-friendly names -----------------------
clean_names <- function(x) {
  x %>%
    as.character() %>%
    gsub(" |/|-", "_", .) %>%
    gsub("[^0-9A-Za-z_]", "", .) %>%
    factor()
}

setwd("J:/MF_Single_Cell_Updated_No_Static")
seur <- readRDS("Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj")

# Clean Seurat metadata --------------------------------------------------
seur@meta.data <- seur@meta.data %>%
  mutate(
    predicted_CellType = clean_names(predicted_CellType)
  )

# Read and clean per-sample fraction data --------------------------------
df_long <- readRDS("celltype.counts.fraction.label.transf.arlotta.5mo.Micro.Bio.only.rds") %>%
  mutate(
    predicted_CellType = factor(clean_names(predicted_CellType)),
    value              = as.numeric(value)
  )

#Build per-cell metadata with Condition ---------------------------------
metadata <- seur@meta.data %>%
  mutate(
    sample.simple      = as.factor(sample.simple),
    Condition = case_when(
      grepl("Microfluidic", sample.simple, ignore.case = TRUE) ~ "Microfluidic",
      grepl("Bioreactor",   sample.simple, ignore.case = TRUE) ~ "Bioreactor",
      TRUE                                                      ~ NA_character_
    )
  ) %>%
  filter(!is.na(Condition)) %>%
  mutate(Condition = factor(Condition, levels = c("Bioreactor", "Microfluidic")))

valid_conditions <- c("Bioreactor", "Microfluidic")

# Loop over each cell type
glmm_results <- map_dfr(
  levels(metadata$predicted_CellType),
  function(ct) {
    df_ct <- metadata %>%
      filter(Condition %in% valid_conditions) %>%
      mutate(
        is_ct = as.integer(predicted_CellType == ct),
        Condition = factor(Condition, levels = valid_conditions)
      )
    
    # Skip if not enough data
    if (n_distinct(df_ct$Condition) < 2 || n_distinct(df_ct$sample) < 2) {
      return(tibble(CellType = ct, p_value = NA_real_))
    }
    
    # Fit full and null mixed logistic regression
    ctl <- glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6),
      check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
      check.conv.grad = .makeCC(action = "ignore", tol = 1e-3)
    )
    
    full_mod <- tryCatch(
      glmer(
        is_ct ~ Condition + (1 | sample),
        data = df_ct,
        family = binomial,
        control = ctl
      ),
      error = function(e) NULL
    )
    
    null_mod <- tryCatch(
      glmer(
        is_ct ~ 1 + (1 | sample),
        data = df_ct,
        family = binomial,
        control = ctl
      ),
      error = function(e) NULL
    )
    
    # Compare with ANOVA if both fit
    pval <- if (!is.null(full_mod) && !is.null(null_mod)) {
      anova(null_mod, full_mod, test = "Chisq")$`Pr(>Chisq)`[2]
    } else {
      NA_real_
    }
    
    tibble(CellType = ct, p_value = pval)
  }
)

# Multiple testing correction (BH)
glmm_results <- glmm_results %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p_value)

glmm_results

# Prepare df_long for plotting ------------------------------------------
df_long <- df_long %>%
  mutate(
    Condition = case_when(
      grepl("Microfluidic", sample, ignore.case = TRUE) ~ "Microfluidic",
      grepl("Bioreactor",   sample, ignore.case = TRUE) ~ "Bioreactor",
      TRUE                                               ~ NA_character_
    )
  ) %>%
  filter(!is.na(Condition)) %>%
  mutate(Condition = factor(Condition, levels = c("Bioreactor", "Microfluidic")))

# Ensure “plots/” directory exists -------------------------------------
plots_dir <- file.path(getwd(), "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Define colors for each cell type -------------------------------------
celltype_colors <- c(
  Astroglia       = "#fb9a99",        
  CFuPN           = "#bebada",           
  CPN             = "#a50f15",             
  Glial_precursors= "#8dd3c7",
  IN_progenitors  = "#fccde5",  
  IP              = "#02818a",             
  Immature_IN     = "#dfc27d",     
  Unspecified_PN  = "#8c510a",  
  aRG             = "#08519C",             
  oRG             = "#f6e8c3",             
  oRG_Astroglia   = "#c7e9c0"   
)

# Loop over cell types to plot & save ----------------------------------
for (ct in levels(df_long$predicted_CellType)) {
  
  ct_data <- df_long %>%
    filter(predicted_CellType == ct)
  
  if (nrow(ct_data) == 0) next
  
  y_ann <- max(ct_data$value, na.rm = TRUE) * 1.1
  
  pval_info <- glmm_results %>%
    filter(CellType == ct) %>%
    pull(adj_p_value)
  
  if (length(pval_info) != 1) pval_info <- NA_real_
  
  signif_label <- case_when(
    pval_info < 0.001 ~ "***",
    pval_info < 0.01  ~ "**",
    pval_info < 0.05  ~ "*",
    TRUE              ~ "ns"
  )
  
  ct_data <- df_long %>%
    filter(predicted_CellType == ct) %>%
    mutate(
      replicate = sample,
      replicate = factor(replicate)
    )
  
  if (nrow(ct_data) == 0) next
  
  ct_data <- ct_data %>%
    arrange(Condition, replicate) %>%
    mutate(replicate = factor(replicate, levels = unique(replicate)))
  
  n_bioreactor <- ct_data %>% filter(Condition == "Bioreactor") %>% pull(replicate) %>% n_distinct()
  vline_pos <- n_bioreactor + 0.5
  
  pval_info <- glmm_results %>%
    filter(CellType == ct) %>%
    pull(adj_p_value)
  if (length(pval_info) != 1) pval_info <- NA_real_
  pval_label <- paste0("Adj. p = ", signif(pval_info, 3))
  y_ann <- max(ct_data$value, na.rm = TRUE) * 1.15
  
  p <- ggplot(ct_data, aes(x = replicate, y = value, fill = Condition)) +
    geom_col(color = "black", width = 0.7, alpha = 0.8) +
    
    scale_fill_manual(values = c(
      Bioreactor   = celltype_colors[[ct]],
      Microfluidic = celltype_colors[[ct]]
    )) +
    
    scale_color_manual(values = c(
      Bioreactor   = "black",
      Microfluidic = "black"
    )) +
    
    geom_vline(
      xintercept = vline_pos,
      linetype   = "dotted",
      color      = "grey40",
      linewidth  = 0.8
    ) +
    annotate(
      "text",
      x = mean(c(1, length(levels(ct_data$replicate)))),
      y = y_ann,
      label = pval_label,
      size = 4,
      fontface = "bold"
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(ct, " Proportion by Condition"),
      x     = "Replicate",
      y     = "Proportion"
    ) +
    theme(
      axis.line       = element_line(linewidth = 0.8, color = "black"),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y     = element_text(size = 10, color = "black"),
      axis.title      = element_text(size = 12, face = "bold"),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  out_file <- file.path(plots_dir, paste0(ct, "_prop_glmm_replicates.pdf"))
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 6, 
    height   = 4,
    units    = "in"
  )
  
  
}

print(glmm_results)
