suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(emmeans)
})

# ============================================================
# Inputs
# ============================================================
rf_all_rds <- "J:/MF and Fetal Data Analysis/project2tuannguy_229RF_MicrofluidicsFetalToQuadrato_RF_AllTri_CC_MT_RP_Excl_3kHVG_Wider_MtryModel_from_fetal_All/quadrato_with_preds.rds"
fetal_rds  <- "J:/MF and Fetal Data Analysis/fetal_dataset_Wang_2025_all.Robj"

out_dir <- "J:/MF and Fetal Data Analysis/FetalTransfer_OppositeYAxis_Props_v5_normalized_labels_recheck"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred_col <- "RF_pred_from_fetal_All"

organoid_condition_col <- "sample.simple"
organoid_sample_col    <- "sample"

fetal_label_col  <- "subclass"
fetal_group_col  <- "Group"
fetal_sample_col <- "donor_id"

# Editable relabel / collapse block
rename_map <- NULL
collapse_map <- list(
  progenitors = c("IPC EN", "Radial glia", "IPC Glia")
)

keep_labels <- NULL
# Example:
# keep_labels <- c("Glutamatergic_neuron", "GABAergic_neuron", "Astrocyte", "progenitors")

organoid_levels <- c("Bioreactor", "Microfluidic")

# Fetal line color only (no shaded band)
fetal_line_color <- "#0072B2"

organoid_colors <- c(
  "Bioreactor"   = "#E69F00",
  "Microfluidic" = "#9E79B9"
)

# ============================================================
# Helpers
# ============================================================
normalize_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[_\\-]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

clean_names <- function(x) {
  x %>%
    as.character() %>%
    gsub(" |/|-", "_", ., perl = TRUE) %>%
    gsub("[^0-9A-Za-z_]", "", ., perl = TRUE) %>%
    factor()
}

recode_labels <- function(x, collapse_map = list(), rename_map = NULL) {
  x <- normalize_label(x)
  
  if (!is.null(rename_map) && length(rename_map) > 0) {
    rename_names <- normalize_label(names(rename_map))
    names(rename_map) <- rename_names
    idx <- x %in% names(rename_map)
    x[idx] <- rename_map[x[idx]]
  }
  
  if (length(collapse_map) > 0) {
    for (new_name in names(collapse_map)) {
      old_vals <- normalize_label(collapse_map[[new_name]])
      x[x %in% old_vals] <- new_name
    }
  }
  
  x
}

sample_totals <- function(df) {
  df %>%
    group_by(sample_id, ConditionGroup) %>%
    summarise(total_cells = n(), .groups = "drop")
}

fit_celltype_glmm <- function(df, ct) {
  all_levels <- c("Bioreactor", "Microfluidic",
                  "First_trimester", "Second_trimester", "Third_trimester")
  
  df_ct <- df %>%
    mutate(
      is_ct = as.integer(CellType == ct),
      ConditionGroup = factor(ConditionGroup, levels = all_levels)
    )
  
  if (n_distinct(df_ct$ConditionGroup) < 2 || n_distinct(df_ct$sample_id) < 2) {
    return(list(
      overall = tibble(CellType = ct, p_value = NA_real_),
      pairwise = tibble()
    ))
  }
  
  ctl <- glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5),
    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
    check.conv.grad = .makeCC(action = "ignore", tol = 1e-3)
  )
  
  full_mod <- tryCatch(
    suppressWarnings(
      glmer(
        is_ct ~ ConditionGroup + (1 | sample_id),
        data = df_ct,
        family = binomial,
        control = ctl
      )
    ),
    error = function(e) NULL
  )
  
  null_mod <- tryCatch(
    suppressWarnings(
      glmer(
        is_ct ~ 1 + (1 | sample_id),
        data = df_ct,
        family = binomial,
        control = ctl
      )
    ),
    error = function(e) NULL
  )
  
  overall_p <- if (!is.null(full_mod) && !is.null(null_mod)) {
    anova(null_mod, full_mod, test = "Chisq")$`Pr(>Chisq)`[2]
  } else {
    NA_real_
  }
  
  if (is.null(full_mod)) {
    return(list(
      overall = tibble(CellType = ct, p_value = overall_p),
      pairwise = tibble()
    ))
  }
  
  emm <- tryCatch(
    emmeans(full_mod, ~ ConditionGroup),
    error = function(e) NULL
  )
  
  if (is.null(emm)) {
    return(list(
      overall = tibble(CellType = ct, p_value = overall_p),
      pairwise = tibble()
    ))
  }
  
  pair_df <- as.data.frame(pairs(emm, adjust = "BH")) %>%
    mutate(
      CellType = ct,
      group1 = vapply(strsplit(as.character(contrast), " - ", fixed = TRUE), `[`, character(1), 1),
      group2 = vapply(strsplit(as.character(contrast), " - ", fixed = TRUE), `[`, character(1), 2),
      label = paste0("Adj. p = ", signif(p.value, 3))
    )
  
  list(
    overall = tibble(CellType = ct, p_value = overall_p),
    pairwise = pair_df
  )
}

# ============================================================
# Load objects
# ============================================================
rf_all <- readRDS(rf_all_rds)
fetal  <- readRDS(fetal_rds)

stopifnot(pred_col %in% colnames(rf_all@meta.data))
stopifnot(organoid_condition_col %in% colnames(rf_all@meta.data))
stopifnot(organoid_sample_col %in% colnames(rf_all@meta.data))

stopifnot(fetal_label_col %in% colnames(fetal@meta.data))
stopifnot(fetal_group_col %in% colnames(fetal@meta.data))
stopifnot(fetal_sample_col %in% colnames(fetal@meta.data))

# ============================================================
# Build organoid metadata
# ============================================================
meta_org <- rf_all@meta.data %>%
  mutate(
    ConditionGroup = case_when(
      .data[[organoid_condition_col]] == "Bioreactor"   ~ "Bioreactor",
      .data[[organoid_condition_col]] == "Microfluidic" ~ "Microfluidic",
      TRUE ~ NA_character_
    ),
    sample_id = as.character(.data[[organoid_sample_col]]),
    CellType_raw = as.character(.data[[pred_col]]),
    source = "organoid"
  ) %>%
  filter(!is.na(ConditionGroup), !is.na(sample_id), !is.na(CellType_raw)) %>%
  mutate(
    CellType_raw = recode_labels(CellType_raw, collapse_map = collapse_map, rename_map = rename_map),
    CellType = clean_names(CellType_raw),
    ConditionGroup = factor(ConditionGroup, levels = organoid_levels)
  )

# ============================================================
# Build fetal metadata
# ============================================================
meta_fetal <- fetal@meta.data %>%
  mutate(
    ConditionGroup = as.character(.data[[fetal_group_col]]),
    sample_id = paste0(.data[[fetal_group_col]], "__", .data[[fetal_sample_col]]),
    CellType_raw = as.character(.data[[fetal_label_col]]),
    source = "fetal"
  ) %>%
  filter(ConditionGroup %in% c("First_trimester", "Second_trimester", "Third_trimester")) %>%
  filter(!is.na(sample_id), !is.na(CellType_raw)) %>%
  mutate(
    CellType_raw = recode_labels(CellType_raw, collapse_map = collapse_map, rename_map = rename_map),
    CellType = clean_names(CellType_raw),
    ConditionGroup = factor(ConditionGroup, levels = c("First_trimester", "Second_trimester", "Third_trimester"))
  )

meta2 <- bind_rows(meta_org, meta_fetal) %>%
  filter(!is.na(ConditionGroup), !is.na(sample_id), !is.na(CellType))

if (!is.null(keep_labels)) {
  meta2 <- meta2 %>% filter(as.character(CellType) %in% keep_labels)
}

# Debug check for collapse
print(table(meta2$CellType))

write.csv(meta2, file.path(out_dir, "filtered_metadata.csv"), row.names = FALSE)

# ============================================================
# Compute percentages per sample and cell type
# ============================================================
sample_props <- meta2 %>%
  group_by(sample_id, ConditionGroup, CellType) %>%
  summarise(ct_cells = n(), .groups = "drop") %>%
  left_join(sample_totals(meta2), by = c("sample_id", "ConditionGroup")) %>%
  mutate(percentage = 100 * ct_cells / total_cells)

write.csv(sample_props, file.path(out_dir, "sample_percentages.csv"), row.names = FALSE)

# ============================================================
# Statistics
# ============================================================
ct_levels <- levels(factor(meta2$CellType))

overall_stats <- map_dfr(ct_levels, function(ct) {
  fit_celltype_glmm(meta2, ct)$overall
})

pairwise_stats <- map_dfr(ct_levels, function(ct) {
  fit_celltype_glmm(meta2, ct)$pairwise
})

write.csv(overall_stats, file.path(out_dir, "overall_glmm_stats.csv"), row.names = FALSE)
write.csv(pairwise_stats, file.path(out_dir, "pairwise_glmm_stats_all.csv"), row.names = FALSE)

bio_mf_stats <- pairwise_stats %>%
  filter(
    (group1 == "Bioreactor" & group2 == "Microfluidic") |
      (group1 == "Microfluidic" & group2 == "Bioreactor")
  )

write.csv(bio_mf_stats, file.path(out_dir, "pairwise_glmm_stats_bio_mf_only.csv"), row.names = FALSE)

# ============================================================
# Plotting
# ============================================================
for (ct in ct_levels) {
  plot_data <- sample_props %>%
    filter(CellType == ct)
  
  organoid_data <- plot_data %>%
    filter(ConditionGroup %in% organoid_levels) %>%
    mutate(sample_id = factor(sample_id, levels = unique(sample_id)))
  
  fetal_summary <- plot_data %>%
    filter(ConditionGroup %in% c("First_trimester", "Second_trimester", "Third_trimester")) %>%
    group_by(ConditionGroup) %>%
    summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      sd_pct = ifelse(n() > 1, sd(percentage, na.rm = TRUE), 0),
      .groups = "drop"
    ) %>%
    mutate(
      ConditionGroup = factor(ConditionGroup, levels = c("First_trimester", "Second_trimester", "Third_trimester"))
    )
  
  if (nrow(organoid_data) == 0 && nrow(fetal_summary) == 0) next
  
  pval_info <- bio_mf_stats %>%
    filter(CellType == ct) %>%
    pull(label)
  
  if (length(pval_info) != 1) pval_info <- NA_character_
  pval_label <- ifelse(is.na(pval_info), "Adj. p = NA", pval_info)
  
  y_max <- max(c(organoid_data$percentage, fetal_summary$mean_pct), na.rm = TRUE)
  y_top <- y_max * 1.20
  
  p <- ggplot() +
    geom_boxplot(
      data = organoid_data,
      aes(x = ConditionGroup, y = percentage, fill = ConditionGroup),
      width = 0.45,
      outlier.shape = NA,
      alpha = 0.75,
      color = "black"
    ) +
    geom_jitter(
      data = organoid_data,
      aes(x = ConditionGroup, y = percentage),
      width = 0.10,
      size = 1.9,
      alpha = 0.7,
      color = "black"
    ) +
    geom_hline(
      data = fetal_summary,
      aes(yintercept = mean_pct),
      linetype = "dotted",
      color = fetal_line_color,
      linewidth = 1
    ) +
    geom_text(
      data = fetal_summary,
      aes(
        x = 2.75,
        y = mean_pct,
        label = ConditionGroup
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.8,
      fontface = "italic",
      color = fetal_line_color
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_top * 0.99,
      label = pval_label,
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = organoid_colors,
      breaks = c("Bioreactor", "Microfluidic"),
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      sec.axis = dup_axis(name = "Percent")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(ct, " proportion by condition"),
      x = NULL,
      y = "Percent"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 90, 10, 10)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(ct, "_organoid_fetal_oppositeYAxis_minusRange.pdf")),
    plot = p,
    width = 6.2,
    height = 4.0,
    units = "in"
  )
}

message("Done. Outputs saved to: ", out_dir)

