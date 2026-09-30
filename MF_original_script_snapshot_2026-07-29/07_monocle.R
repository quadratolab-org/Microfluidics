
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("batchelor", "BiocSingular"))
BiocManager::install(c("DelayedArray", "SingleCellExperiment", "SummarizedExperiment"))



# Load devtools
library(devtools)

# Install monocle3 from Cole Trapnell's GitHub
devtools::install_github("cole-trapnell-lab/monocle3")
devtools::install_github("satijalab/seurat-wrappers")

library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)


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
  "oRG/Astroglia" = "#c7e9c0",
  "Unknown" = '#d9d9d9'
)


setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS("Final_joined_obj_label_transfer_Arlotta.Robj")
seur <- subset(seur, subset = data.set == "quadrato" )
head(seur)
seur

Idents(seur)<- 'predicted_CellType'
unique(Idents(seur))
seur <- subset(seur, subset = predicted_CellType != "Unknown")


# Subset to equal number of cells per condition
set.seed(123)

conditions <- c("Bioreactor", "Microfluidic", "Static")

# Find the minimum number of cells in any sample
sample_counts <- table(seur$sample)  
min_cells <- min(sample_counts)

# Downsample each sample to the smallest group size
subset_seur <- subset(seur, cells = unlist(lapply(names(sample_counts), function(s) {
  sample(Cells(seur)[seur$sample == s], min_cells)
})))

# Check the new sample distribution
table(subset_seur$sample)

setwd("J:/MF_Single_Cell_Updated_No_Static")
seur <- readRDS('Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj')

# STEP 2: Convert to Monocle3 CellDataSet
cds <- as.cell_data_set(seur)

cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Transfer cluster and UMAP embeddings (if not already done)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

colnames(colData(cds))

# STEP 3: Plot cells to help identify where to root the trajectory
# Generate the plot
p <- plot_cells(
  cds,
  color_cells_by = "predicted_CellType",
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 0.8,  # Reduce the point size
  show_trajectory_graph = TRUE
) +
  scale_color_manual(values = celltype_colors) +
  theme_minimal(base_size = 14) +  # Seurat-like minimal theme
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),  # Adjust legend size
    plot.title = element_text(hjust = 0.5, size = 16),  # Title in the center
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)  # Clean margin for a more polished look
  ) +
  theme(
    # Adjust label appearance
    strip.text = element_text(size = 16, face = "bold", color = "black"),  # Increase label size and make it bold
      )
print(p)

# Save the plot as a PDF
setwd("J:/MF_Single_Cell_Updated_No_Static/Monocle")
ggsave("monocle_trajectory_umap.pdf", plot = p, width = 8, height = 6)  # Customize size as needed


# ⚠️ STEP 4: Manually pick a root cell 
progenitor_cells <- rownames(seur@meta.data %>% 
                               filter(predicted_CellType == "aRG"))

# Use one of these cells as the root (can use plot_cells() to refine this)
root_cell <- progenitor_cells[1]

# STEP 5: Order cells in pseudotime from that root
cds <- order_cells(cds, root_cells = root_cell)

# STEP 6: Visualize pseudotime
p2<-plot_cells(cds, color_cells_by = "pseudotime")
colnames(colData(cds))
p2
ggsave("monocle_trajectory_pseudotime_umap.pdf", plot = p2, width = 8, height = 6)  # Customize size as needed


# STEP 7: Kolmogorov–Smirnov test to compare pseudotime across conditions

pseudotime_df <- data.frame(
  pseudotime = pseudotime(cds),
  condition = colData(cds)$sample.simple
) %>% filter(!is.na(pseudotime))

# Perform pairwise one-sided KS tests
condition_pairs <- combn(conditions, 2, simplify = FALSE)
ks_test_results <- list()

for (pair in condition_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  test_result <- ks.test(
    pseudotime_df$pseudotime[pseudotime_df$condition == group1],
    pseudotime_df$pseudotime[pseudotime_df$condition == group2],
    alternative = "less"  # Test if group1 is accelerated compared to group2
  )
  
  ks_test_results[[paste(group1, "vs", group2)]] <- test_result
}

# STEP 8: Print KS test results
ks_test_results



# Remove non-finite (NA or Inf) values from the pseudotime data
pseudotime_df_clean <- pseudotime_df[!is.na(pseudotime_df$pseudotime) & is.finite(pseudotime_df$pseudotime), ]

# Create a density plot with improved aesthetics
p3<-ggplot(pseudotime_df_clean, aes(x = pseudotime, fill = condition, color = condition)) +
  geom_density(alpha = 0.6, linewidth = 1.2) +  # Use linewidth instead of size
  labs(
    title = "Pseudotime Distribution Across Conditions",
    x = "Pseudotime", 
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  scale_fill_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  scale_color_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",  # Place the legend at the top
    legend.box = "horizontal",
    panel.grid.major = element_blank(),   # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on the left and right
  scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on the top and bottom
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1)  # Add a vertical dashed line for reference

p3
ggsave("monocle_density_plots.pdf", plot = p3, width = 8, height = 6)  # Customize size as needed


# Create a density plot with improved aesthetics and separate panels for each condition
ggplot(pseudotime_df_clean, aes(x = pseudotime, fill = condition, color = condition)) +
  geom_density(alpha = 0.6, linewidth = 1.2) +  # Use linewidth instead of size
  labs(
    title = "Pseudotime Distribution Across Conditions",
    x = "Pseudotime", 
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  scale_fill_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  scale_color_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",  # Place the legend at the top
    legend.box = "horizontal",
    panel.grid.major = element_blank(),   # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on the left and right
  scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on the top and bottom
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # Add a vertical dashed line for reference
  facet_wrap(~ condition, scales = "free_y")  # Separate plots for each condition




######################################################## test these out and figure out what I was trying to do hahaha 04/22/25

#Extract UMAP coordinates and metadata
umap_df <- as.data.frame(reducedDims(cds)$UMAP)

# Rename the unnamed columns (usually V1, V2) to UMAP_1 and UMAP_2
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Add metadata (like predicted_CellType)
umap_df$celltype <- colData(cds)$celltype
umap_df$sample <- colData(cds)$sample
# Optional: Clean factor levels for better labeling
umap_df$celltype <- factor(umap_df$celltype)

# Calculate centroids for labeling
library(dplyr)
label_positions <- umap_df %>%
  group_by(celltype) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# ✨ Prettier UMAP Plot
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text(data = label_positions, aes(label = celltype), 
            color = "black", size = 4.5, fontface = "bold", check_overlap = TRUE) +
  theme_minimal(base_size = 15) +
  labs(
    title = "UMAP Colored by Predicted Cell Types",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_viridis_d(option = "plasma")




#######################################################################################################

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = predicted_CellType)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text(
    data = label_positions,
    aes(label = predicted_CellType),
    color = "black",
    size = 4.5,
    fontface = "bold",
    check_overlap = TRUE
  ) +
  theme_minimal(base_size = 15) +
  labs(
    title = "UMAP Colored by Predicted Cell Types",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_manual(values = celltype_colors)




############################## split by samples

label_positions <- umap_df %>%
  group_by(celltype, sample) %>%
  summarize(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text(
    data = label_positions,
    aes(x = UMAP_1, y = UMAP_2, label = celltype),
    color = "black",
    size = 4.5,
    fontface = "bold",
    check_overlap = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~sample) +
  theme_minimal(base_size = 15) +
  labs(
    title = "UMAP by Sample (Colored by Predicted Cell Type)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_color_manual(values = celltype_colors)



####################################################################################################
####################################################################################################
####################################################################################################
##################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#02818a','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#dd3497','#67000d','#fa9fb5','#d9d9d9')

# running monocle on trevino pcw 24 dataset plus ours

library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("N:/JoJo/Team Microfluidic-JP/seurat.object/Trevino.PCW24.static.micro.bio")
load("object.neurons.rg.ipc.cca.Robj")
setwd("~/Microfluidics McCain Collab/single_cell_run/Monocle")
object.neuron.rg.ipc.lognorm<- object.cca
head(object.neuron.rg.ipc.lognorm)
rm(object.cca)
unique(object.neuron.rg.ipc.lognorm$class)
unique(object.neuron.rg.ipc.lognorm$sample.celltype)
unique(object.neuron.rg.ipc.lognorm$simple.class)

object.neuron.rg.ipc.lognorm$sample <- sub("-.*", "",
                                           object.neuron.rg.ipc.lognorm$sample.celltype
)

object.neuron.rg.ipc.lognorm$celltype <- sub(".*-", "", 
                                             object.neuron.rg.ipc.lognorm$sample.celltype
)
unique(object.neuron.rg.ipc.lognorm$celltype)
unique(object.neuron.rg.ipc.lognorm$sample)

object.neuron.rg.ipc.lognorm <- JoinLayers(object.neuron.rg.ipc.lognorm)
object.neuron.rg.ipc.lognorm

# Subset to equal number of cells per condition
set.seed(123)
conditions <- c("bio", "micro", "static", 'fetal')

# Identify the sample ID column 
sample_counts <- table(object.neuron.rg.ipc.lognorm$sample)  

# Find the minimum number of cells in any sample
min_cells <- min(sample_counts)

# Downsample each sample to the smallest group size
subset_seur <- subset(object.neuron.rg.ipc.lognorm, cells = unlist(lapply(names(sample_counts), function(s) {
  sample(Cells(object.neuron.rg.ipc.lognorm)[object.neuron.rg.ipc.lognorm$sample == s], min_cells)
})))
subset_seur
# Check the new sample distribution
table(subset_seur$sample)

DimPlot(subset_seur, reduction = 'umap.cca', group.by = 'celltype', split.by = 'sample', cols = cols)

# Convert to Monocle3 CellDataSet
cds <- as.cell_data_set(subset_seur)

cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition = FALSE) #troubleshooting line
cds <- learn_graph(cds)

colnames(colData(cds))

#troubleshooting code
cds <- preprocess_cds(cds, num_dim = 50)   #troubleshooting code
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.metric = "cosine", umap.min_dist = 0.1, umap.n_neighbors = 30)   #troubleshooting code
cds <- cluster_cells(cds, resolution = 1e-2)  #troubleshooting code
cds <- learn_graph(cds, use_partition = FALSE) #troubleshooting code
############### end ###########


library(RColorBrewer)
celltype_colors <- brewer.pal(n = 12, name = "Set3")[1:11]

p <- plot_cells(
  cds,
  color_cells_by = "celltype",
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 0.8,  
  show_trajectory_graph = TRUE
) +
  scale_color_manual(values = celltype_colors) +
  theme_minimal(base_size = 14) +  # Seurat-like minimal theme
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),  # Adjust legend size
    plot.title = element_text(hjust = 0.5, size = 16),  # Title in the center
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)  
  ) +
  theme(
    strip.text = element_text(size = 16, face = "bold", color = "black"),  
  )

print(p)
ggsave("fetal_monocle_trajectory_umap.pdf", plot = p, width = 8, height = 6)  





# ⚠️ Manually pick a root cell 
progenitor_cells <- rownames(subset_seur@meta.data %>% 
                               filter(celltype == "Early_RG"))

# Use one of these cells as the root (can use plot_cells() to refine this)
root_cell <- progenitor_cells[1]

# Order cells in pseudotime from that root
cds <- order_cells(cds, root_cells = root_cell)

# Visualize pseudotime
plot_cells(cds, color_cells_by = "pseudotime")
colnames(colData(cds))



# Kolmogorov–Smirnov test to compare pseudotime across conditions
pseudotime_df <- data.frame(
  pseudotime = pseudotime(cds),
  condition = colData(cds)$sample
) %>% filter(!is.na(pseudotime))

# Perform pairwise one-sided KS tests
condition_pairs <- combn(conditions, 2, simplify = FALSE)
ks_test_results <- list()

for (pair in condition_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  test_result <- ks.test(
    pseudotime_df$pseudotime[pseudotime_df$condition == group1],
    pseudotime_df$pseudotime[pseudotime_df$condition == group2],
    alternative = "less"  # Test if group1 is accelerated compared to group2
  )
  
  ks_test_results[[paste(group1, "vs", group2)]] <- test_result
}


ks_test_results



# Remove non-finite (NA or Inf) values from the pseudotime data
pseudotime_df_clean <- pseudotime_df[!is.na(pseudotime_df$pseudotime) & is.finite(pseudotime_df$pseudotime), ]

conditions <- c("bio", "micro", "static", 'fetal')


# Create a density plot with improved aesthetics
ggplot(pseudotime_df_clean, aes(x = pseudotime, fill = condition, color = condition)) +
  geom_density(alpha = 0.6, linewidth = 1.2) +  # Use linewidth instead of size
  labs(
    title = "Pseudotime Distribution Across Conditions",
    x = "Pseudotime", 
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  scale_fill_manual(values = c("static" = "#1b9e77", "micro" = "#d95f02", "bio" = "#7570b3", 'fetal'= 'grey')) +
  scale_color_manual(values = c("static" = "#1b9e77", "micro" = "#d95f02", "bio" = "#7570b3", 'fetal'= 'grey')) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",  # Place the legend at the top
    legend.box = "horizontal",
    panel.grid.major = element_blank(),   # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on the left and right
  scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on the top and bottom
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1)  # Add a vertical dashed line for reference


# Create a density plot with improved aesthetics and separate panels for each condition
ggplot(pseudotime_df_clean, aes(x = pseudotime, fill = condition, color = condition)) +
  geom_density(alpha = 0.6, linewidth = 1.2) +  # Use linewidth instead of size
  labs(
    title = "Pseudotime Distribution Across Conditions",
    x = "Pseudotime", 
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  scale_fill_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  scale_color_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",  # Place the legend at the top
    legend.box = "horizontal",
    panel.grid.major = element_blank(),   # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on the left and right
  scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on the top and bottom
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # Add a vertical dashed line for reference
  facet_wrap(~ condition, scales = "free_y")  # Separate plots for each condition













































