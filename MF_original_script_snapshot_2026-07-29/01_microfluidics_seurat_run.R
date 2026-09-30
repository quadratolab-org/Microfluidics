############
#Microfluidics single cell
###########

library(BiocManager)
library(glmGamPoi)
library(dplyr)
#library(Seurat, lib.loc = 'C:/CustomR')
library(Seurat)
library(viridis)
library(ggplot2)
library(stringr)
library(reshape2)
library(DESeq2)

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

seur3 <- subset(seur3, subset = sample.simple != "Static")


Idents(seur3)<- 'predicted_CellType'
unique(Idents(seur3))
seur4 <- subset(seur3, subset = predicted_CellType != "Unknown")
seur4
#seur4
#An object of class Seurat 
#44206 features across 38373 samples within 1 assay 
#Active assay: RNA (44206 features, 2000 variable features)
#3 layers present: data, counts, scale.data
#3 dimensional reductions calculated: pca, integrated.cca, umap.cca


min_cells <- min(table(seur4$sample))
set.seed(123)


cells_to_keep <- seur4@meta.data %>%
  mutate(cell_id = Cells(seur4)) %>%  
  group_by(sample) %>%
  sample_n(min_cells) %>%
  pull(cell_id)


seur4_downsampled <- subset(seur4, cells = cells_to_keep)

seur4_downsampled
table(seur4_downsampled$sample)

saveRDS(seur4_downsampled, file= 'Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj')


# previous object from Jo
#load('final_cca_object.Robj')
#head(object.cca)
#object.cca
#head(object.cca)


object.cca.meta <- object.cca@meta.data
head(object.cca.meta)
class(object.cca.meta)

object.cca.meta<- object.cca.meta %>%
  mutate(across('sample', str_replace_all, '_1', ''))

object.cca.meta<- object.cca.meta %>%
  mutate(across('sample', str_replace_all, '_2', ''))

object.cca.meta<- object.cca.meta %>%
  mutate(across('sample', str_replace_all, '_3', ''))

object.cca.meta<- object.cca.meta %>%
  mutate(across('sample', str_replace_all, '_4', ''))

colnames(object.cca.meta)[which(names(object.cca.meta) == "sample")] <- "sample_simple"


object.cca <- AddMetaData(object.cca, metadata = object.cca.meta, col.name = 'sample_simple')

object.cca.meta <- object.cca@meta.data

#what name are in sample column
categories <- unique(object.cca.meta$sample_simple) 
categories
numberOfCategories <- length(categories)
numberOfCategories



load('object.join.final2.Robj')


head(object.join)
object.join.tail<- tail(object.join)
class(object.join.tail)
write.csv(object.join.tail, 'object.join.tail.csv')


head(joined.object)
joined.object

joined.meta<- joined.object@meta.data
head(joined.meta)
class(joined.meta)
write.table(joined.meta, file = 'joined.meta.txt')

#converting object.join which is log normalized as a scanpy anndata readable object
library(Seurat)
library(SeuratData)
library(SeuratDisk)

SaveH5Seurat(joined.object, filename = "joined.object.h5Seurat")
Convert("joined.object.h5Seurat", dest = "h5ad")

joined.object[["RNA3"]] <- as(object = joined.object[["RNA"]], Class = "Assay")
DefaultAssay(joined.object) <- "RNA3"
joined.object[["RNA"]] <- NULL
joined.object <- RenameAssays(object = joined.object, RNA3 = 'RNA')

# load('object.join.final2.Robj') ### joined layer object 090524

setwd('C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run/seurat figures')


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


pdf(file = 'split_test.pdf', width=48, height=10, useDingbats=FALSE)
DimPlot(object.join, group.by = 'celltype1', pt.size = 1, cols = cols)
DimPlot(object.join, group.by = 'celltype1', split.by = 'sample', pt.size = 1, cols = cols)
DimPlot(object.join, group.by = 'celltype1', split.by = 'sample.simple', pt.size = 1, cols = cols)

dev.off()

features= c('PGK1', 'ARCN1', 'GORASP1')
FeaturePlot(object.join, features = c('PGK1', 'ARCN1', 'GORASP1'))
FeaturePlot(object.join, features = c('SATB2', 'NRP1', 'CAV1'))

features2= c('SATB2', 'BCL11B', 'NEUROD6', 'CAMK2A', 'SYN1', 'GRIA2', 'SLC17A7')


RidgePlot(object.join, features = features, ncol = 2)

pdf(file = 'Bhaduri_stress_markers_VlnPlot.pdf', width=48, height=10, useDingbats=FALSE)
VlnPlot(object.join, features = features2, group.by = 'sample.simple', alpha = 0, cols = cols )
dev.off()

pdf(file = 'DotPlot_neuronal_genes.pdf', width=24, height=10, useDingbats=FALSE)
DotPlot(object.join, features = features2, group.by = "sample.simple") + RotatedAxis()
dev.off()

getwd()

VlnPlot(object.join, features = features, split.by = 'celltype1', idents = c('CFuPN', 'CPN', 'aRG'), group.by = 'sample.simple', alpha = 0, cols = cols )

head(object.join)

png(file = 'CCA_Merged_Cell_Types.png', width=800, height=800)
DimPlot(object.cca, group.by = 'celltype1', split.by = 'data.set', pt.size = 1, cols = cols)
dev.off()

#FeaturePlot(object.cca, features = c('FOXG1', 'FEZF2', 'SATB2', 'HOPX', 'SOX2', 'PAX6', 'TOP2A', 'TTR',
                                # 'NEUROD6', 'EOMES', 'GAD1', 'BCL11B',
                                # 'GAD67', 'RORB', 'AUTS2'))

FeaturePlot(object.cca, features = c('ROBO1', 'ROBO2', 'ROBO3', 'ROBO4', 'PLXNA2', 'PLXNA4', 'SATB2', 'FEZF2', 'BCL11B',
                                     'NRP1', 'NRP2', 'L1CAM', 'DCC, UNC5' ,'DOCK3', 'NRCAM'))

features2= c('SATB2', 'BCL11B', 'NEUROD6', 'CAMK2A', 'SYN1', 'GRIA2', 'SLC17A7')

object.join.ExN<- subset(object.join, idents= c('CFuPN', 'CPN'))
DotPlot(object.join.ExN, features = features2, group.by = "sample.simple") + RotatedAxis()


#################################################################################

########### Label Transfer from Uzquiano 5 month old cortical Organoid ##########

###############################################################################

# Ensure 'data.set' is a factor
object.join$data.set <- as.factor(object.join$data.set)

# Split into reference (Arlotta) and query (Quadrato)
ref_obj <- subset(object.join, subset = data.set == "arlotta")
query_obj <- subset(object.join, subset = data.set == "quadrato")

# Find transfer anchors
anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30)

# Transfer cell type labels
query_obj <- TransferData(anchorset = anchors, refdata = ref_obj$CellType, dims = 1:30)
head(query_obj) 
class(query_obj)

# Merge back into the original Seurat object
object.join$predicted_CellType <- NA  # Create an empty column
object.join$predicted_CellType[rownames(query_obj)] <- query_obj$predicted.id
# View results
head(object.join@meta.data)

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



pdf(file = 'Label_transfer_umap.pdf', width=48, height=10, useDingbats=FALSE)
DimPlot(seur, group.by = 'predicted_CellType', split.by = 'sample.simple', pt.size = 1, cols = celltype_colors)
dev.off()

pdf(file = 'Label_transfer_umap_sample_split.pdf', width=48, height=10, useDingbats=FALSE)
DimPlot(object.join, group.by = 'predicted_CellType', split.by = 'sample', pt.size = 1, cols = celltype_colors)
dev.off()

pdf(file = 'IN_Feature_PLOTLabel_transfer_umap.pdf', width=48, height=36, useDingbats=FALSE)
FeaturePlot(seur, features = c('HOPX', 'RPS6', 'MOXD1', 'GFAP'), split.by = 'sample.simple', pt.size = 1, label.size = 6 )
dev.off()

par(mfrow = c(1, 2))

DimPlot(seur, group.by = 'predicted_CellType', split.by = 'sample.simple', pt.size = 1, cols = celltype_colors)
FeaturePlot(seur, features = c('EOMES', 'NEUROD6', 'SOX2', 'DLX5'), split.by = 'sample.simple', pt.size = 1, label.size = 6 )



head(object.join)
Idents(object.join)<- 'celltype1'
counts <- unlist(lapply(levels(Idents(object.join)), function(ident)
{length(WhichCells(object.join, idents = ident))}))
names(counts) <- levels(Idents(object.join))
counts

total_cells <- sum(counts)
print(total_cells)

#Choroid_Plexus                                CFuPN                            Astroglia                                  aRG 
#1268                                 4306                                 1353                                 2702 



#Cycling_oRG/All        Cycling_Immature_IN/oRG/Glial_Progen            Cycling_IPC/oRG              Mix_IPC/oRG/Immature_IN 
#681                                 3090                                 1884                                 1684 



#Mix_oRG/oRG-Astro/Astroglia          immature_IN                        oRG-Astroglia                       Unspecified_PN 
#2123                                 3326                                10474                                11693 



#IPC                                  CPN 
#5644                                20618 

head(object.join)
Idents(object.join)<- 'predicted_CellType'
counts2 <- unlist(lapply(levels(Idents(object.join)), function(ident)
{length(WhichCells(object.join, idents = ident))}))
names(counts2) <- levels(Idents(object.join))
counts2

total_cells2 <- sum(counts2)
print(total_cells2)

#bebada           #08519C           #02818a           #f6e8c3         #a50f15       #d9d9d9       #c7e9c0              #fb9a99       #8c510a           #dfc27d
#CFuPN              aRG               IP              oRG              CPN          Unknown    oRG/Astroglia        Astroglia   Unspecified PN      Immature IN 
#3168             3577             3542             2296            14995              686             4673              732             4110              839 

#fccde5             #8dd3c7
#IN progenitors    Glial precursors 
#232                  209 

object.join <- subset(object.join, subset = data.set == "quadrato")
tail(object.join)
object.join
Assays(object.join)

saveRDS(object.join, file= 'Final_joined_obj_label_transfer_Arlotta.Robj')
seur<-readRDS('Final_joined_obj_label_transfer_Arlotta.Robj')
head(seur)
seur
rm(anchors)

################################################################
############# barplots ##################################
######################################################
library(reshape2)
library(tidyr)

counts = as.matrix(table(seur$predicted_CellType, seur$sample)) #set so that first column is clusters/celltypes, second is organoids, or however you want to split up the x-axis
counts = t(t(counts)/colSums(counts)) #transform from raw counts to percentages, to normalize for library size

counts<- as.data.frame(counts)

# Rename columns for clarity
colnames(counts) <- c("predicted_CellType", "sample", "value")

head(counts)
class(counts)

counts <- counts %>%
  mutate(Condition = case_when(
    grepl("Static", sample) ~ "Static",
    grepl("Microfluidic", sample) ~ "Microfluidic",
    grepl("Bioreactor", sample) ~ "Bioreactor"
  ))

write.csv(counts,"celltype.counts.fraction.label.transf.arlotta.5mo.csv", row.names = FALSE)
saveRDS(counts, file = "celltype.counts.fraction.label.transf.arlotta.5mo.rds")

#### AMI, mutual information analysis from Uzquiano
library(aricode)
library(ggplot2)
#AMI: cite Vinh et al "Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance"

seur.meta <- seur@meta.data
head(seur.meta)
static <- seur.meta[seur.meta$sample.simple == 'Static', ]
bioreactor <- seur.meta[seur.meta$sample.simple == 'Bioreactor', ]
microfluidic <- seur.meta[seur.meta$sample.simple == 'Microfluidic', ]

#Mutual information between clusters and organoids
adjusted_mi_static = AMI(static$predicted_CellType, static$sample)
adjusted_mi_bioreactor = AMI(bioreactor$predicted_CellType, bioreactor$sample)
adjusted_mi_microfluidic = AMI(microfluidic$predicted_CellType, microfluidic$sample)


# Enhanced bar plot with spacing between conditions
p1 <- ggplot(counts, aes(fill = predicted_CellType, y = value, x = sample)) + 
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors) +  # Custom colors
  facet_grid(~Condition, scales = "free_x", space = "free") +  # Facets by condition with free space
  labs(
    title = "Proportion Cell Types Across Samples",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text
    axis.title.x = element_text(size = 12, face = "bold"),  # Bold x-axis label
    axis.title.y = element_text(size = 12, face = "bold"),  # Bold y-axis label
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "bottom",  # Position legend at the bottom
    legend.title = element_text(size = 10, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 9)  # Adjust legend text size
  ) +
  # Add annotations for AMI scores
  annotate("text", x = 1, y = 1.05, label = paste("AMI:", round(adjusted_mi_static, 5)), size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = 5, y = 1.05, label = paste("AMI:", round(adjusted_mi_microfluidic, 5)), size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = 9, y = 1.05, label = paste("AMI:", round(adjusted_mi_bioreactor, 5)), size = 4, fontface = "bold", hjust = 0)


# Show the plot
p1
getwd()
setwd("~/Microfluidics McCain Collab/single_cell_run/seurat figures/Module_Scores_Arlotta_Transfer")

ggsave("Proportion_of_Cell_Types_Across_Samples .tiff", plot = last_plot(), width = 10, height = 8, dpi = 300, device = "tiff")



################## Pairwise comparisons between cell types across conditions


# Load necessary libraries
library(tidyverse)
library(reshape2)
library(rstatix)  # for pairwise comparison with stats

# Load the data
df_long <- readRDS("celltype.counts.fraction.label.transf.arlotta.5mo.rds")
# Make sure numeric column is numeric
df_long$value <- as.numeric(df_long$value)

# Get all unique cell types
cell_types <- unique(df_long$predicted_CellType)

# Pairwise t-tests per predicted_CellType
pairwise_results <- map_df(cell_types, function(ct) {
  
  ct_data <- df_long %>% filter(predicted_CellType == ct)
  
  ct_data <- ct_data %>%
    mutate(Condition = sub("_[0-9]+$", "", sample))
  
  pairwise <- pairwise.wilcox.test(ct_data$value, ct_data$Condition, p.adjust.method = "BH")
  
  broom::tidy(pairwise) %>%
    rename(adj.p.value = p.value) %>%
    mutate(predicted_CellType = ct)
})

# Preview results
head(pairwise_results)


library(tidyverse)
library(ggpubr)

# Make sure consistent condition labels
condition_levels <- c("Static", "Bioreactor", "Microfluidic")

# Recreate output folder if needed
dir.create("plots", showWarnings = FALSE)

library(ggpubr)  # Needed for stat_pvalue_manual

for (ct in unique(df_long$predicted_CellType)) {
  
  # Clean and prep data
  ct_data <- df_long %>%
    filter(predicted_CellType == ct) %>%
    mutate(
      Condition = sub("_[0-9]+$", "", sample),
      Condition = factor(Condition, levels = condition_levels)
    )
  
  # Prepare comparisons
  ct_comparisons <- pairwise_results %>%
    filter(predicted_CellType == ct) %>%
    mutate(
      group1 = factor(group1, levels = levels(ct_data$Condition)),
      group2 = factor(group2, levels = levels(ct_data$Condition)),
      y.position = max(ct_data$value, na.rm = TRUE) * 1.05 + row_number() * 0.02
    )
  
  ct_comparisons$Condition <- factor(levels(ct_data$Condition)[1], levels = levels(ct_data$Condition))
  
  # Reorder the factor levels
  ct_data$Condition <- factor(ct_data$Condition, levels = c("Bioreactor", "Microfluidic", "Static"))
  
  p <- ggplot(ct_data, aes(x = Condition, y = value)) +
    # Bar plot with single cell type color
    stat_summary(
      fun = mean,
      geom = "bar",
      fill = celltype_colors[ct],
      color = "black",
      width = 0.6,
      size = 0.8,
      alpha = 0.8
    ) +
    
    # Error bars
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      size = 0.6
    ) +
    
    # Jitter points
    geom_jitter(color = 'black', width = 0.2, size = 3, alpha = 0.7) +
    
    theme_classic(base_size = 10) +
    theme(
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    
    labs(
      title = paste(ct, "Across Conditions"),
      x = NULL,
      y = "Proportion"
    ) +
    
    stat_compare_means(
      comparisons = list(
        c("Bioreactor", "Microfluidic"), 
        c("Bioreactor", "Static"),
        c("Microfluidic", "Static")
      ),
      method = "t.test",
      label = "p.signif",
      size = 4,
      hide.ns = TRUE
    ) +
    
    stat_compare_means(
      comparisons = list(
        c("Bioreactor", "Microfluidic"), 
        c("Bioreactor", "Static"),
        c("Microfluidic", "Static")
      ),
      method = "t.test",
      label = "p.format",
      size = 3,
      vjust = 2
    )
  
  
  
  
  
  # Save
  ggsave(filename = paste0("plots/", ct, "_barplot.pdf"), plot = p, width = 4.5, height = 4, units = 'in', dpi = 300)
}




############################### redoing cell type porportion stats using mixed-effect logistic regression model as in Paulsen et al.,
# using the seur that has unknwon removed; generated the counts table from the previous bar plot attempt
# You can treat condition (ie. BR, MF,ST) as a categorical fixed effect and include it directly in the model

install.packages("lme4")
install.packages("dplyr")
install.packages("broom.mixed")

library(lme4)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(ggpubr)

# Assume `metadata` has:
# - predicted_CellType
# - sample.simple
# - Condition (e.g., Static, Bioreactor, Microfluidic)

# Ensure factors
seur$sample <- as.factor(seur$sample)
seur$sample.simple <- as.factor(seur$sample.simple)
seur$predicted_CellType <- as.factor(seur$predicted_CellType)
metadata <- seur@meta.data

# Create condition column
metadata$Condition <- case_when(
  grepl("Static", metadata$sample) ~ "Static",
  grepl("Microfluidic", metadata$sample) ~ "Microfluidic",
  grepl("Bioreactor", metadata$sample) ~ "Bioreactor"
)
metadata$Condition <- factor(metadata$Condition, levels = c("Bioreactor", "Microfluidic", "Static"))
head(metadata)
# Prepare pairwise comparisons
comparisons <- combn(levels(metadata$Condition), 2, simplify = FALSE)

# Store pairwise GLMM results
pairwise_glmm_results <- list()

for (ct in levels(seur$predicted_CellType)) {
  ct_df <- metadata
  ct_df$celltype_binary <- as.integer(ct_df$predicted_CellType == ct)
  
  for (comp in comparisons) {
    sub_df <- ct_df %>% filter(Condition %in% comp)
    sub_df$Condition <- droplevels(sub_df$Condition)
    
    model_full <- glmer(celltype_binary ~ Condition + (1 | sample), 
                        data = sub_df, 
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))
    
    model_null <- glmer(celltype_binary ~ (1 | sample), 
                        data = sub_df, 
                        family = binomial,
                        control = glmerControl(optimizer = "bobyqa"))
    
    model_comp <- anova(model_null, model_full)
    
    pairwise_glmm_results[[paste(ct, comp[1], comp[2], sep = "_")]] <- data.frame(
      CellType = ct,
      group1 = comp[1],
      group2 = comp[2],
      p_value = model_comp$`Pr(>Chisq)`[2]
    )
  }
}

######################################################################################################
############ this one was working!!! #######################################

glmm_results <- list()

for (ct in levels(seur$predicted_CellType)) {
  
  # Create binary outcome for this cell type
  metadata$celltype_binary <- as.integer(metadata$predicted_CellType == ct)
  
  # Full model with Condition
  model_full <- glmer(celltype_binary ~ Condition + (1 | sample.simple), 
                      data = metadata, 
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))
  
  # Reduced model (no Condition)
  model_null <- glmer(celltype_binary ~ (1 | sample.simple), 
                      data = metadata, 
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))
  
  # Compare models
  model_comp <- anova(model_null, model_full)
  
  # Store results
  glmm_results[[ct]] <- data.frame(
    CellType = ct,
    p_value = model_comp$`Pr(>Chisq)`[2]
  )
}

# Combine and adjust p-values
glmm_df <- do.call(rbind, glmm_results)
glmm_df$adj_p_value <- p.adjust(glmm_df$p_value, method = "BH")
glmm_df <- glmm_df %>% arrange(adj_p_value)

print(glmm_df)


##################################################################################################

# Combine and adjust
pairwise_df <- bind_rows(pairwise_glmm_results)
pairwise_df$adj_p_value <- p.adjust(pairwise_df$p_value, method = "BH")
head(metadata)
# Compute per-sample proportions
cell_counts <- metadata %>%
  group_by(sample, Granular.Cells) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cell_counts <- cell_counts %>%
  mutate(Condition = case_when(
    grepl("Static", sample) ~ "Static",
    grepl("Microfluidic", sample) ~ "Microfluidic",
    grepl("Bioreactor", sample) ~ "Bioreactor"
  )) %>%
  mutate(Condition = factor(Condition, levels = c("Bioreactor", "Microfluidic", "Static")))
head(cell_counts)
# Plot
for (ct in unique(cell_counts$Granular.Cells)) {
  ct_data <- cell_counts %>% filter(Granular.Cells == ct)
  ymax <- max(ct_data$prop, na.rm = TRUE)
  
  # Get pairwise p-values for current cell type
  signif_data <- pairwise_df %>%
    filter(CellType == ct) %>%
    mutate(
      y.position = ymax + 0.05 + row_number() * 0.03,
      significance = case_when(
        adj_p_value < 0.001 ~ "***",
        adj_p_value < 0.01 ~ "**",
        adj_p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  p <- ggplot(ct_data, aes(x = Condition, y = prop)) +
    stat_summary(fun = mean, geom = "bar", fill = c( "#7570b3", "#d95f02","#1b9e77"), color = "black", width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.6, color = "black") +
    
    stat_pvalue_manual(
      signif_data,
      label = "significance",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 4
    ) +
    
    theme_classic(base_size = 13) +
    labs(
      title = paste0(ct, " Cell Proportion"),
      x = NULL,
      y = "Proportion"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold")
    )
  
  # Save
  ggsave(filename = paste0("plots/", ct, "_glmm_pairwise_barplot.pdf"), plot = p, width = 5, height = 4)
}


getwd()













library(dplyr)
library(lme4)
library(ggplot2)
library(ggpubr)
library(purrr)


seur$sample <- as.factor(seur$sample)
seur$sample.simple <- as.factor(seur$sample.simple)
seur$predicted_CellType <- as.factor(seur$predicted_CellType)

metadata <- seur@meta.data

# Create condition variable
metadata$Condition <- case_when(
  grepl("Static", metadata$sample) ~ "Static",
  grepl("Microfluidic", metadata$sample) ~ "Microfluidic",
  grepl("Bioreactor", metadata$sample) ~ "Bioreactor"
)
metadata$Condition <- factor(metadata$Condition, levels = c("Bioreactor", "Microfluidic", "Static"))






# List of pairwise condition comparisons
comparisons <- combn(levels(metadata$Condition), 2, simplify = FALSE)

# Store results
glmm_all_results <- list()
pairwise_glmm_results <- list()

for (ct in levels(metadata$predicted_CellType)) {
  ct_df <- metadata
  ct_df$celltype_binary <- as.integer(ct_df$predicted_CellType == ct)
  
  # 1. Full model (all 3 conditions)
  model_full <- glmer(celltype_binary ~ Condition + (1 | sample), 
                      data = ct_df, 
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))
  
  model_null <- glmer(celltype_binary ~ (1 | sample), 
                      data = ct_df, 
                      family = binomial,
                      control = glmerControl(optimizer = "bobyqa"))
  
  model_comp <- anova(model_null, model_full)
  
  glmm_all_results[[ct]] <- data.frame(
    CellType = ct,
    p_value = model_comp$`Pr(>Chisq)`[2]
  )
  
  # 2. Pairwise comparisons
  for (comp in comparisons) {
    sub_df <- ct_df %>% filter(Condition %in% comp)
    sub_df$Condition <- droplevels(sub_df$Condition)
    
    pw_full <- glmer(celltype_binary ~ Condition + (1 | sample), 
                     data = sub_df, 
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))
    
    pw_null <- glmer(celltype_binary ~ (1 | sample), 
                     data = sub_df, 
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))
    
    pw_comp <- anova(pw_null, pw_full)
    
    pairwise_glmm_results[[paste(ct, comp[1], comp[2], sep = "_")]] <- data.frame(
      CellType = ct,
      group1 = comp[1],
      group2 = comp[2],
      p_value = pw_comp$`Pr(>Chisq)`[2]
    )
  }
}







# Initialize result lists
glmm_all_results <- list()
pairwise_glmm_results <- list()

# Define pairwise combinations of conditions
comparisons <- combn(levels(metadata$Condition), 2, simplify = FALSE)

# Loop over cell types
for (ct in levels(metadata$predicted_CellType)) {
  
  ct_df <- metadata
  ct_df$celltype_binary <- as.integer(ct_df$predicted_CellType == ct)
  
  # --- Full model with all 3 conditions ---
  try({
    model_full <- glmer(
      celltype_binary ~ Condition + (1 | sample), 
      data = ct_df, 
      family = binomial,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    
    model_null <- glmer(
      celltype_binary ~ (1 | sample), 
      data = ct_df, 
      family = binomial,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    
    # Check if model is singular
    if (isSingular(model_full)) {
      message(paste("Singular fit detected in full model for:", ct, "- refitting with fixed effects only."))
      
      model_full <- glm(celltype_binary ~ Condition, data = ct_df, family = binomial)
      model_null <- glm(celltype_binary ~ 1, data = ct_df, family = binomial)
    }
    
    model_comp <- anova(model_null, model_full, test = "Chisq")
    
    glmm_all_results[[ct]] <- data.frame(
      CellType = ct,
      p_value = model_comp$`Pr(>Chi)`[2]
    )
  }, silent = TRUE)
  
  # --- Pairwise GLMMs ---
  for (comp in comparisons) {
    sub_df <- ct_df %>% filter(Condition %in% comp)
    sub_df$Condition <- droplevels(sub_df$Condition)
    
    try({
      pw_full <- glmer(
        celltype_binary ~ Condition + (1 | sample), 
        data = sub_df, 
        family = binomial,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      
      pw_null <- glmer(
        celltype_binary ~ (1 | sample), 
        data = sub_df, 
        family = binomial,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      
      # Check for singular fit
      if (isSingular(pw_full)) {
        message(paste("Singular fit in pairwise for:", ct, paste(comp, collapse = "_"), "- switching to fixed-effects model."))
        
        pw_full <- glm(celltype_binary ~ Condition, data = sub_df, family = binomial)
        pw_null <- glm(celltype_binary ~ 1, data = sub_df, family = binomial)
      }
      
      pw_comp <- anova(pw_null, pw_full, test = "Chisq")
      
      pairwise_glmm_results[[paste(ct, comp[1], comp[2], sep = "_")]] <- data.frame(
        CellType = ct,
        group1 = comp[1],
        group2 = comp[2],
        p_value = pw_comp$`Pr(>Chi)`[2]
      )
    }, silent = TRUE)
  }
}



















###############################################################
  
########## Adding Module Scores for various MSigDB gene sets #########################
      library(data.table)
      
      seur    
      setwd("~/Microfluidics McCain Collab/single_cell_run")

     
      # Read the TSV file
      msigdb_glycolysis <- fread("HALLMARK_GLYCOLYSIS.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_glycolysis[msigdb_glycolysis[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_glycolysis <- unlist(strsplit(gene_symbols_row[, 2], ","))      
      
      
      
      
      # Read the TSV file
      msigdb_mtor <- fread("HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_mtor[msigdb_mtor[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_mtor <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
      
      
      # Read the TSV file
      msigdb_ox_phos <- fread("HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_ox_phos[msigdb_ox_phos[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_ox_phos <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
    
      # Read the TSV file
      msigdb_hypox <- fread("HALLMARK_HYPOXIA.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_hypox[msigdb_hypox[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_hypox <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
      
      
      # Read the TSV file
      msigdb_ros <- fread("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_ros[msigdb_ros[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_ros <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
      
      
      # Read the TSV file
      msigdb_apoptosis <- fread("HALLMARK_APOPTOSIS.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_apoptosis[msigdb_apoptosis[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_apoptosis <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
      
      # Read the TSV file
      msigdb_estrogen <- fread("HALLMARK_ESTROGEN_RESPONSE_EARLY.v2024.1.Hs.tsv", header = FALSE, sep = "\t", data.table = FALSE)
      
      # Extract the row where the first column is "GENE_SYMBOLS"
      gene_symbols_row <- msigdb_estrogen[msigdb_estrogen[, 1] == "GENE_SYMBOLS", ]
      
      # Convert the second column (column index 2) into a vector by splitting it on commas
      msigdb_estrogen <- unlist(strsplit(gene_symbols_row[, 2], ","))  
      
      
   gene_sets<- list(msigdb_apoptosis, msigdb_estrogen, msigdb_glycolysis, msigdb_hypox, msigdb_mtor, msigdb_ox_phos, msigdb_ros)
   names(gene_sets) <- c("Apoptosis_Score", "Estrogen_Score", "Glycolysis_Score", "Hypoxia_Score", "mTOR_Score", "Oxidative_Phosphorylation_Score", "ROS_Score")   
   
   
seur<- AddModuleScore(seur, features = gene_sets, name = 'ModuleScore' )

head(seur)      
colnames(seur@meta.data)[grep("ModuleScore", colnames(seur@meta.data))] <- names(gene_sets)
      
      
      
  
      
      #################################################      
      #################  #########################
      #####################################################       

# HNOCA meta data and glycolysis and other scores


setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")
      
#loading in hnoca meta with glycolysis scores as well

hnoca_meta <- read.csv("obs.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))

head(hnoca_meta)
tail(hnoca_meta)
class(hnoca_meta)

rownames(hnoca_meta) <- hnoca_meta$X  # Set row names
hnoca_meta <- hnoca_meta[, -which(names(hnoca_meta) == "X")]  # Remove the original column


library(dplyr)

hnoca_meta_subset <- hnoca_meta %>% filter(type == "query")

#object.join.hnoca <- AddMetaData(object.join, metadata = hnoca_meta)
head(hnoca_meta_subset)
colnames(hnoca_meta_subset)

write.csv(hnoca_meta_subset, file = 'hnoca_meta_subset_quadrato.csv')
hnoca_meta_subset<-read.csv('hnoca_meta_subset_quadrato.csv', header = T) ### changed column "annot_level_4_rev2_prediction_knn_filtered_by_uncert.0.5" to 'Cell_ID' manually as column OG heading was difficult to wrangel in R
########### script to try for Glyzolysis analysis using hnoca labels#######################

glycolysis_avg <- hnoca_meta_subset %>%
  group_by(batch, Cell_ID) %>%  
  summarise(Average_Glycolysis = mean(Hallmark_Glycolysis, na.rm = TRUE), .groups = "drop")


head(glycolysis_avg)
class(glycolysis_avg)
write.csv(glycolysis_avg, "glycolysis_avg.csv", row.names = F)

unique(glycolysis_avg$Cell_ID)
glycolysis_avg_neurons <- glycolysis_avg %>% filter(Cell_ID == "Dorsal Telencephalic Neuron NT-VGLUT") 
head(glycolysis_avg_neurons)


glycolysis_avg_neurons <- glycolysis_avg_neurons %>% arrange(batch)


df2 <- glycolysis_avg_neurons %>%
  mutate(Batch_Group = case_when(
    grepl("Static", batch) ~ "Static",
    grepl("Microfluidic", batch) ~ "Microfluidic",
    grepl("Bioreactor", batch) ~ "Bioreactor"
  ))

head(df2)

anova_result <- aov(Average_Glycolysis ~ Batch_Group, data = df2)
summary(anova_result)

TukeyHSD(anova_result)


################################################################################################################
library(ggplot2)
library(ggpubr)
library(dplyr)


ggplot(df2, aes(x = Batch_Group, y = Average_Glycolysis, fill = Batch_Group)) +
  
  # Boxplot: Filled with the same colors as dots, black border
  geom_boxplot(width = 0.5, color = "black", size = 1, alpha = 0.6, outlier.shape = NA) +
  
  # Jittered scatter points (larger size)
  geom_jitter(aes(color = Batch_Group), width = 0.2, size = 4, alpha = 0.8) +
  
  # Custom colors for both box fill & dots
  scale_fill_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  scale_color_manual(values = c("Static" = "#1b9e77", "Microfluidic" = "#d95f02", "Bioreactor" = "#7570b3")) +
  
  # GraphPad Prism-style theme
  theme_classic() +
  theme(
    axis.line = element_line(size = 1, color = "black"),  # Thicker axis lines
    axis.text = element_text(size = 14, color = "black"), # Larger axis labels
    axis.title = element_text(size = 16, face = "bold"),  # Bold axis titles
    legend.position = "none",  # Hide legend
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Centered title
  ) +
  
  # Labels
  labs(
    title = "Glycolysis Score Across Conditions",
    x = "Batch Group",
    y = "Average Glycolysis Score in VGLUT Neurons"
  ) +
  
  # Statistical comparisons (p-values & asterisks)
  stat_compare_means(comparisons = list(
    c("Static", "Microfluidic"), 
    c("Static", "Bioreactor"),
    c("Microfluidic", "Bioreactor")
  ), method = "t.test", label = "p.signif", size = 6, hide.ns = TRUE) +  # Asterisks for significance
  
  # Add numerical p-values below asterisks
  stat_compare_means(comparisons = list(
    c("Static", "Microfluidic"), 
    c("Static", "Bioreactor"),
    c("Microfluidic", "Bioreactor")
  ), method = "t.test", label = "p.format", size = 4, vjust = 2)  # Numerical p-values slightly below


getwd()
ggsave("Glycolysis_Score_Neurons.pdf", plot = last_plot(), width = 10, height = 8, dpi = 300)






#################################################################################

##################       pseudo bulk     #################################

#################################################################


# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(DESeq2)
library(tibble)
# Assume `seurat_obj` is your Seurat object with metadata containing 'Condition' and 'Sample_ID'
head(object.join)

object.join <- subset(object.join, subset = data.set == "quadrato" )
# Step 1: Aggregate counts (Pseudobulk) ---------------------------------

# Extract sparse matrix safely and reshape for pseudobulk


library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)

# Ensure metadata has rownames as a column
object.join@meta.data <- object.join@meta.data %>%
  rownames_to_column(var = "Cell_ID")

# Extract the counts matrix as a sparse matrix
counts_matrix <- GetAssayData(object.join, assay = "RNA", slot = "counts")  

# Convert sparse matrix to a long-format dataframe
counts_long <- as.data.frame(summary(counts_matrix))  
colnames(counts_long) <- c("Gene_Index", "Cell_Index", "Count")  

# Convert indices to actual gene and cell names
counts_long$Gene <- rownames(counts_matrix)[counts_long$Gene_Index]
counts_long$Cell_ID <- colnames(counts_matrix)[counts_long$Cell_Index]

# Join metadata
pseudobulk_counts <- counts_long %>%
  left_join(object.join@meta.data %>% select(Cell_ID, sample.simple, sample), by = "Cell_ID") %>%
  group_by(Gene, sample, sample.simple) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = Count, values_fill = 0)

# **Make gene names unique**
pseudobulk_counts$Gene <- make.unique(pseudobulk_counts$Gene)

# Set unique gene names as row names
pseudobulk_counts <- pseudobulk_counts %>%
  column_to_rownames(var = "Gene")


# Step 2: Create DESeq2 dataset ---------------------------------

# Check column names of pseudobulk_counts
print(colnames(pseudobulk_counts))

# Check row names of col_data
print(rownames(col_data))

# Ensure col_data only contains matching samples from pseudobulk_counts
col_data <- col_data[rownames(col_data) %in% colnames(pseudobulk_counts), , drop = FALSE]

# Ensure pseudobulk_counts only contains matching samples from col_data
pseudobulk_counts <- pseudobulk_counts[, colnames(pseudobulk_counts) %in% rownames(col_data)]
col_data <- col_data[colnames(pseudobulk_counts), , drop = FALSE]
# Check if they match now
print(all(colnames(pseudobulk_counts) == rownames(col_data))) # Should return TRUE

# Convert to DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts, colData = col_data, design = ~ sample.simple)

# Step 3: Run Differential Expression Analysis ------------------
dds <- estimateSizeFactors(dds, type = "poscounts")  # Use alternative size factor estimation
dds <- DESeq(dds)  # Proceed with DESeq2 analysis

# Step 4: Extract Results ---------------------------------------
colData(dds)
res.BR.MF <- results(dds, contrast = c("sample.simple", "Bioreactor", "Microfluidic"))
res.BR.MF <- res.BR.MF[order(res.BR.MF$pvalue), ]  # Order by significance

# View top differentially expressed genes
head(res.BR.MF)




library(ggplot2)


# Perform PCA on variance-stabilized transformed counts
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = "sample.simple", returnData = TRUE)  # Get PCA data

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = sample.simple)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples",
       x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% Variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% Variance")) +
  theme(text = element_text(size = 14))




####### labeling top genes ##########################


library(ggplot2)
library(dplyr)
library(ggrepel)  # Load the package

# Convert results to a data frame
res_df <- as.data.frame(res.BR.MF)
res_df$Gene <- rownames(res_df)

# Add significance categories
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Select the top 10 most significant genes based on padj
top_genes <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(10)  # Select top 10 genes

# Plot Volcano with labels
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(text = element_text(size = 14)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # Log2FC thresholds
  geom_text_repel(data = top_genes, aes(label = Gene), size = 5, max.overlaps = 10, box.padding = 0.3, segment.size = 0.3)  # Label top genes



########################################################
########## repeating for MF v Static ###########
##########################################################

res.MF.ST <- results(dds, contrast = c("sample.simple", "Microfluidic", "Static"))
res.MF.ST <- res.MF.ST[order(res.MF.ST$pvalue), ]  # Order by significance

# View top differentially expressed genes
head(res.MF.ST)

# Perform PCA on variance-stabilized transformed counts
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = "sample.simple", returnData = TRUE)  # Get PCA data

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = sample.simple)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples",
       x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% Variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% Variance")) +
  theme(text = element_text(size = 14))

####### labeling top genes ##########################

library(ggplot2)
library(dplyr)
library(ggrepel)  # Load the package

# Convert results to a data frame
res_df <- as.data.frame(res.MF.ST)
res_df$Gene <- rownames(res_df)

# Add significance categories
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Select the top 10 most significant genes based on padj
top_genes <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(10)  # Select top 10 genes

dim(res_df)
head(res_df)
# Plot Volcano with labels
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(text = element_text(size = 14)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # Log2FC thresholds
  geom_text_repel(data = top_genes, aes(label = Gene), size = 5, max.overlaps = 10, box.padding = 0.3, segment.size = 0.3)  # Label top genes



# Identify mitochondrial genes
mito_genes <- rownames(dds)[grepl("^MT-", rownames(dds), ignore.case = TRUE)]

# Calculate percentage of reads mapping to mitochondrial genes per sample
mito_fraction <- colSums(counts(dds)[mito_genes, ]) / colSums(counts(dds))

# Plot mitochondrial content across samples
library(ggplot2)
ggplot(data.frame(Sample = colnames(dds), MitoFraction = mito_fraction),
       aes(x = Sample, y = MitoFraction)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Mitochondrial Gene Fraction Per Sample", y = "Fraction of Reads")




############################### SFARI Analysis #############################################

setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")

#loading in SFARI data to see which DEGs are expressed across conditions

sfari <- read.csv("SFARI-Gene_genes_01-13-2025release_03-11-2025export.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
head(sfari)
dim(sfari)

##########################################################################
#################### Satterstorm ASD DEG #################################
##########################################################################
setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")
asd <- read.csv("ASD.102.genes.only.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
class(asd)
str(asd)
head(seur)
seur_CPN <- subset(seur, subset = predicted_CellType %in% c( "CPN"))
seur_CFuPN <- subset(seur, subset = predicted_CellType %in% c("CFuPN"))
head(seur_CPN)



#
################# Still need to troubleshoot pseudobulk 03212025!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##################################
####


# Extract metadata and ensure conditions are properly labeled
seur_CPN$Condition <- factor(seur_CPN$sample.simple, levels = c("Bioreactor", "Microfluidic", "Static"))
seur_CPN$sample <- as.factor(seur_CPN$sample)

# Create pseudobulk matrix by aggregating counts per sample and condition
pseudobulk_counts <- AverageExpression(
  seur_CPN,
  assays = "RNA",  # Adjust if using SCT or another assay
  features = NULL,  # Default is all features
  group.by = "sample"
)$RNA

# Add metadata for samples
sample_metadata <- seur_CPN@meta.data %>%
  distinct(sample, Condition, predicted_CellType) %>%
  arrange(sample, predicted_CellType)

rownames(sample_metadata) <- sample_metadata$sample
rownames(sample_metadata) <- gsub("_", "-", rownames(sample_metadata))

head(sample_metadata)
head(pseudobulk_counts)

any(!apply(pseudobulk_counts, 1:2, is.integer))
any(pseudobulk_counts < 0)
str(sample_metadata)
str(pseudobulk_counts)

# Convert pseudobulk_counts to integers
pseudobulk_counts_int <- round(as.matrix(pseudobulk_counts))

# Ensure column names of count matrix match row names of sample metadata
colnames(pseudobulk_counts_int) <- gsub("-", "_", colnames(pseudobulk_counts_int))
rownames(sample_metadata) <- gsub("-", "_", sample_metadata$sample)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts_int,
  colData = sample_metadata,
  design = ~ Condition
)

# Perform DESeq normalization and analysis
dds <- DESeq(dds)

# Bioreactor vs Static
res_bioreactor_vs_static <- results(dds, contrast = c("Condition", "Bioreactor", "Static"))
# Microfluidic vs Static
res_microfluidic_vs_static <- results(dds, contrast = c("Condition", "Microfluidic", "Static"))
# Bioreactor vs Microfluidic
res_bioreactor_vs_microfluidic <- results(dds, contrast = c("Condition", "Bioreactor", "Microfluidic"))

# Filter significant DEGs
deg_bioreactor_vs_static <- as.data.frame(res_bioreactor_vs_static) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

deg_microfluidic_vs_static <- as.data.frame(res_microfluidic_vs_static) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

deg_bioreactor_vs_microfluidic <- as.data.frame(res_bioreactor_vs_microfluidic) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# Save results to CSV
setwd("~/Microfluidics McCain Collab/single_cell_run/seurat figures/Module_Scores_Arlotta_Transfer")

write.csv(deg_bioreactor_vs_static, "DEGs_Bioreactor_vs_Static.csv", row.names = TRUE)
write.csv(deg_microfluidic_vs_static, "DEGs_Microfluidic_vs_Static.csv", row.names = TRUE)
write.csv(deg_bioreactor_vs_microfluidic, "DEGs_Bioreactor_vs_Microfluidic.csv", row.names = TRUE)


# Example: Volcano plot for Bioreactor vs Static
volcano_data <- as.data.frame(res_bioreactor_vs_static)
volcano_data$Significant <- volcano_data$padj < 0.05

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Bioreactor vs Static", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  scale_color_manual(values = c("gray", "red"))



##################################################################

########################## Using MAST to calc DEGs for now 03 21 2025 #####################################

#################################################################



if (!requireNamespace("MAST", quietly = TRUE)) {
  install.packages("BiocManager")  # Install BiocManager if not already installed
  BiocManager::install("MAST")     # Install MAST
}
library(MAST)
library(ggrepel)

head(seur)
unique(Idents(seur))

seur_BR_MF <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Microfluidic'))
seur_BR_MF_exc <- subset(seur_BR_MF, subset = predicted_CellType %in% c('CFuPN', 'CPN', 'aRG', 'IP', 'oRG', 'oRG/Astroglia', 'Unspecified PN'))
seur_BR_MF
seur_BR_MF_exc
head(seur_BR_MF)
Idents(seur_BR_MF_exc) <- "sample.simple"

# Perform DEG analysis between conditions
degs <- FindMarkers(
  object = seur_BR_MF_exc,
  ident.1 = "Bioreactor",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA", 'predicted_CellType'),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs <- degs %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Select top genes for labeling (e.g., top 10 up and top 10 down)
top_genes <- degs %>%
  filter(significance != "Not Significant") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 20)

# Volcano plot
library(ggrepel)

volcano_plot <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +  # Scatter plot with transparency
  scale_color_manual(
    values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey80")
  ) +  # Custom color palette
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +  # P-value threshold line
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.8) +  # Fold-change thresholds
    
  labs(
    title = "Volcano Plot: Differentially Expressed Genes Bioreactor vs Microfluidic",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centered and bold title
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  xlim(-2.5, 2.5) +  # Adjust based on your data
  ylim(0, 50)    # Adjust based on your data

# Print plot
print(volcano_plot)

# Save as PDF
ggsave("Volcano_Plot_DEGs_BR_MF_exc.pdf", plot = volcano_plot, width = 10, height = 8, dpi = 300)




list_genes<- as.list(asd)
head(list_genes)
names(list_genes) <- c("X", "genes")
# Convert rownames of DEGs to a column for merging
degs$gene <- rownames(degs)

# Cross the list with DEGs
crossed_genes <- merge(
  x = degs,
  y = list_genes,
  by.x = "gene",   # Column in DEGs to match
  by.y = "genes"       # Column in the provided list to match
)

# View the crossed genes
print(crossed_genes)
rownames(crossed_genes) <- crossed_genes$gene

print(degs)
dim(degs)

# Subset the data for crossed genes
crossed_genes <- degs[rownames(degs) %in% list_genes$genes, ]

volcano_plot <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  # Make crossed genes larger and more visible
  geom_point(
    aes(color = ifelse(rownames(degs) %in% list_genes$genes, "Crossed", "Other")),
    alpha = ifelse(rownames(degs) %in% list_genes$genes, 1, 0.4),  # Brighter for crossed genes
    size = ifelse(rownames(degs) %in% list_genes$genes, 3, 1.5)  # Larger for crossed genes
  ) +
  # Set custom colors
  scale_color_manual(values = c("Crossed" = "red", "Other" = "grey")) +
    # Custom color palette
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +  # P-value threshold line
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.8) +  # Fold-change thresholds
  # Highlight and label crossed genes
  geom_text_repel(
    data = crossed_genes[
      abs(crossed_genes$avg_log2FC) > 0.25 & 
        -log10(crossed_genes$p_val_adj) > -log10(0.05),
    ],
    aes(label = rownames(crossed_genes)[
      abs(crossed_genes$avg_log2FC) > 0.25 & 
        -log10(crossed_genes$p_val_adj) > -log10(0.05)
    ]),,
    size = 4,
    color = "darkred",  # Custom text color
    box.padding = 0.6,
    point.padding = 0.5,
    segment.color = "black",
    max.overlaps = 10
  ) +
  
  # Add plot labels and titles
  labs(
    title = "Volcano Plot of ASD Genes",
    x = "Average Log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Gene Category"
  ) +
  
  # Minimal theme with customizations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # Remove grid lines
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )+
  ylim(0, 50)+
xlim(-2.5, 2.5)

print(volcano_plot)
getwd()
# Save the plot as a PDF
ggsave("Crossed_Genes_Volcano_Plot_Exc.pdf", plot = volcano_plot, width = 10, height = 8, dpi = 300)


##########################################################################################################################################
########################## Run it back, but downsample so we have equal number of celltypes for the MAST analysis ########################
##########################################################################################################################################


head(seur)
DimPlot(seur, group.by = 'predicted_CellType', cols = celltype_colors)
table(seur@meta.data$predicted_CellType)


# Group cells by cell type
cell_groups <- table(seur@meta.data$predicted_CellType)

# Find the minimum number of cells among the available cell types
min_cells <- min(cell_groups)

# Downsample cells for each cell type 
set.seed(123)  # For reproducibility
downsampled_cells <- unlist(lapply(names(cell_groups), function(ct) {
  # Check if the cell type exists in the data
  if (ct %in% seur@meta.data$predicted_CellType) {
    # Subset cells of the current cell type
    cells_of_type <- WhichCells(seur, ident = ct)
    
    # Check if there are enough cells to sample
    if (length(cells_of_type) >= min_cells) {
      sample(cells_of_type, size = min_cells, replace = FALSE)
    } else {
      cells_of_type  # Use all available cells if fewer than min_cells
    }
  }
}))

# Subset Seurat object to the downsampled cells
seur_downsampled <- subset(seur, cells = downsampled_cells)

# Confirm the new Seurat object has equal representation of cell types (for available ones)
table(seur_downsampled@meta.data$predicted_CellType)


head(seur_downsampled)




##########################################################################################################################################
########################## Run it back, but for only neurons for the MAST analysis ########################
##########################################################################################################################################

FeaturePlot(seur, features = c('OSTN', 'GREM2', 'BMF', 'LVRN', 'LAMB4'), split.by = 'sample.simple')


library(MAST)
library(ggrepel)

seur <- readRDS('Final_seur_quadrato_removed_unknown.Robj')

seur
head(seur)
seur_BR_MF <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Microfluidic'))
seur_BR_MF_neurons <- subset(seur_BR_MF, subset = predicted_CellType %in% c( "CFuPN", 'CPN'))

seur_BR_MF
seur_BR_MF_neurons

head(seur_BR_MF_neurons)
Idents(seur_BR_MF_neurons) <- "sample.simple"

# Perform DEG analysis between conditions
degs <- FindMarkers(
  object = seur_BR_MF_neurons,
  ident.1 = "Bioreactor",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs <- degs %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs <- degs %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

output_dir <- "MAST_DEGs"
dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs, "neuron_only_MAST_DEGs_BR_MF.csv", row.names = TRUE)

# Select top genes for labeling (e.g., top 10 up and top 10 down)
top_genes <- degs %>%
  filter(significance != "Not Significant") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 20)

# Volcano plot
library(ggrepel)

volcano_plot <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +  # Scatter plot with transparency
  scale_color_manual(
    values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey80")
  ) +  # Custom color palette
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +  # P-value threshold line
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.8) +  # Fold-change thresholds
  labs(
    title = "Volcano Plot: Differentially Expressed Genes Bioreactor vs Microfluidic Neurons",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centered and bold title
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  xlim(-2.5, 2.5) +  # Adjust based on your data
  ylim(0, 50)    # Adjust based on your data


# Print plot
print(volcano_plot)

# Save as PDF
setwd("~/Microfluidics McCain Collab/single_cell_run/seurat figures")
ggsave("Volcano_Plot_DEGs_BR_MF_neurons_only.pdf", plot = volcano_plot, width = 10, height = 8, dpi = 300)

#######################################################################################
####################### crossing with Satterstorm ASD genes; 102 list ####################
#####################################################################################

list_genes<- as.list(asd)
head(list_genes)
names(list_genes) <- c("X", "genes")
# Convert rownames of DEGs to a column for merging
degs$gene <- rownames(degs)

# Cross the list with DEGs
crossed_genes <- merge(
  x = degs,
  y = list_genes,
  by.x = "gene",   # Column in DEGs to match
  by.y = "genes"       # Column in the provided list to match
)

# View the crossed genes
print(crossed_genes)
rownames(crossed_genes) <- crossed_genes$gene
class(crossed_genes)

# Sort by log fold change (descending order)
crossed_genes <- crossed_genes[order(-crossed_genes$avg_log2FC), ]

# Save the table as a CSV file
write.csv(crossed_genes, file = "crossed_asd_genes_sorted_BR_vs_MF.csv", row.names = TRUE)


# Subset the data for crossed genes
crossed_genes <- degs[rownames(degs) %in% list_genes$genes, ]

volcano_plot2 <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  # Highlight crossed genes prominently
  geom_point(
    aes(color = ifelse(rownames(degs) %in% list_genes$genes, "Crossed", "Other")),
    alpha = 0.8,
    size = ifelse(rownames(degs) %in% list_genes$genes, 3.5, 2)  # Larger for crossed genes
  ) +
  # Custom color palette
  scale_color_manual(
    values = c("Crossed" = "#E41A1C", "Other" = "grey80"),
    labels = c("Crossed" = "Highlighted Genes", "Other" = "All Other Genes")
  ) +
  # Add dashed threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.8) +
  # Add gene names for crossed genes
  geom_text_repel(
    data = degs[rownames(degs) %in% list_genes$genes, ],  # Subset to crossed genes
    aes(label = rownames(crossed_genes)),  # Gene names
    size = 4, 
    color = "black",  # Label color
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = 15  # Adjust to avoid excessive overlap
  ) +
  # Add plot labels and titles
  labs(
    title = "Volcano Plot: Differentially Expressed Genes",
    subtitle = "Bioreactor vs Microfluidic Neurons",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Category"
  ) +
  # Minimal theme with clean customizations
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()
  ) +
  xlim(-2.5, 2.5) +  # Adjust range as needed
  ylim(0, 50)        # Maintain the same y-axis limit for consistency

print(volcano_plot2)
setwd("~/Microfluidics McCain Collab/single_cell_run/seurat figures")
ggsave("Crossed_Genes_Volcano_Plot_Neurons_only.pdf", plot = volcano_plot2, width = 10, height = 8, dpi = 300)





#########################################################################
#crossing with E/L RGs as in Revah et al., 2022

###########################################################################################################################
custom_gene_sets <- list(
  "Human_LRG" = human_LRG,
  "Mouse_ERG" = mouse_erg,
  "Mouse_LRG" = mouse_lrg
)

head(custom_gene_sets)
# Convert rownames of DEGs to a column for merging
degs$gene <- rownames(degs)

# Cross the list with DEGs
crossed_genes <- merge(
  x = degs,
  y = custom_gene_sets,
  by.x = "gene",   # Column in DEGs to match
  by.y = "Gene Name"       # Column in the provided list to match
)

# View the crossed genes
print(crossed_genes)
rownames(crossed_genes) <- crossed_genes$gene
class(crossed_genes)

# Sort by log fold change (descending order)
crossed_genes <- crossed_genes[order(-crossed_genes$avg_log2FC), ]

# Save the table as a CSV file
write.csv(crossed_genes, file = "crossed_asd_genes_sorted_BR_vs_MF.csv", row.names = TRUE)


# Subset the data for crossed genes
crossed_genes <- degs[rownames(degs) %in% list_genes$genes, ]

volcano_plot <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  # Highlight crossed genes prominently
  geom_point(
    aes(color = ifelse(rownames(degs) %in% list_genes$genes, "Crossed", "Other")),
    alpha = 0.8,
    size = ifelse(rownames(degs) %in% list_genes$genes, 3.5, 2)  # Larger for crossed genes
  ) +
  # Custom color palette
  scale_color_manual(
    values = c("Crossed" = "#E41A1C", "Other" = "grey80"),
    labels = c("Crossed" = "Highlighted Genes", "Other" = "All Other Genes")
  ) +
  # Add dashed threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.8) +
  # Add gene names for crossed genes
  geom_text_repel(
    data = degs[rownames(degs) %in% list_genes$genes, ],  # Subset to crossed genes
    aes(label = rownames(crossed_genes)),  # Gene names
    size = 4, 
    color = "black",  # Label color
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = 15  # Adjust to avoid excessive overlap
  ) +
  # Add plot labels and titles
  labs(
    title = "Volcano Plot: Differentially Expressed Genes",
    subtitle = "Bioreactor vs Microfluidic Neurons",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Category"
  ) +
  # Minimal theme with clean customizations
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()
  ) +
  xlim(-2.5, 2.5) +  # Adjust range as needed
  ylim(0, 50)        # Maintain the same y-axis limit for consistency

print(volcano_plot)

ggsave("Crossed_Genes_Volcano_Plot_Neurons_only.pdf", plot = volcano_plot, width = 10, height = 8, dpi = 300)

getwd()




############### Gene Ontology using DEGs from MAST analysis; doing for neurons only and then for all cells; progenitors only? IPC only?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)

# Read in your DEGs CSV — assuming genes as rownames
degs <- read.csv("neuron_only_MAST_DEGs_BR_MF.csv", row.names = 1)

head(degs)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs <- degs %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs$ENTREZID, degs$significance)

# Run compareCluster for GO:BP enrichment
ego_compare <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot results
dotplot(ego_compare, showCategory=15, font.size=10) +
  ggtitle("GO Biological Process Enrichment (Up/Downregulated DEGs)")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results.csv", row.names = FALSE)

pdf("output/GO_BP_dotplot.pdf", width=10, height=8)
print(dotplot(ego_compare, showCategory=15, font.size=10) +
        ggtitle("GO Biological Process Enrichment (Up/Downregulated DEGs)"))
dev.off()

ggsave("output/GO_BP_dotplot.png", width=10, height=8)







