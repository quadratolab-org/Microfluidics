
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
library(MAST)
library(ggrepel)


# using the downsampled, unknown removed object 042525 


setwd("~/Microfluidics McCain Collab/single_cell_run")
#seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur <- readRDS('Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj')

rm(seur2, seur3, seur4, seur4_downsampled)
seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
#[1] 13
length(rb_genes)  # how many ribosomal genes removed
#[1] 108
length(genes_to_remove)  # total genes removed
#[1] 121
dim(seur)  # new dimension of Seurat object
# [1] 44085 21160

seur_BR_MF <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Microfluidic'))
seur_BR_MF_neurons <- subset(seur_BR_MF, subset = predicted_CellType %in% c( "CFuPN", 'CPN'))


Idents(seur_BR_MF_neurons) <- "sample.simple"

# Perform DEG analysis between conditions
degs <- FindMarkers(
  object = seur_BR_MF_neurons,
  ident.1 = "Bioreactor",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size
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

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs, "neuron_only_MAST_DEGs_BR_MF.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=10, font.size=10) +
  ggtitle("GO Biological Process Enrichment Bioreactor vs. Microfluidic Neurons Only")

# 723 upregulated; 1260 downregulated BR vs MF

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 15 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nBioreactor vs. Microfluidic Neurons Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top15_UpDown_BRvsMF.pdf", p, width = 8, height = 7)


#########################################################################################################################
############################################################################
###################################################################################################################

# Static vs. Microfluidic Neurons

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown.Robj')

seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712

seur_ST_MF <- subset(seur, subset = sample.simple %in% c( "Static", 'Microfluidic'))
seur_ST_MF_neurons <- subset(seur_ST_MF, subset = predicted_CellType %in% c( "CFuPN", 'CPN'))


Idents(seur_ST_MF_neurons) <- "sample.simple"

# Perform DEG analysis between conditions
degs2 <- FindMarkers(
  object = seur_ST_MF_neurons,
  ident.1 = "Static",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs2 <- degs2 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs2 <- degs2 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs2, "neuron_only_MAST_DEGs_ST_MF.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs2 <- read.csv("neuron_only_MAST_DEGs_ST_MF.csv", row.names = 1)

head(degs2)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs2), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs2 <- degs2 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs2$ENTREZID, degs2$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=15, font.size=10) +
  ggtitle("GO Biological Process Enrichment Static vs. Microfluidic Neurons Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_STvsMF.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nStatic vs. Microfluidic Neurons Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_STvsMF.pdf", p, width = 8, height = 7)



##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

#### Bioreactor vs Static Neurons


# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712

seur_BR_ST <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Static'))
seur_BR_sT_neurons <- subset(seur_BR_ST, subset = predicted_CellType %in% c( "CFuPN", 'CPN'))


Idents(seur_BR_sT_neurons) <- "sample.simple"

# Perform DEG analysis between conditions
degs3 <- FindMarkers(
  object = seur_BR_sT_neurons,
  ident.1 = "Bioreactor",          
  ident.2 = "Static",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs3 <- degs3 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs3 <- degs3 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs3, "neuron_only_MAST_DEGs_BR_ST.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs3 <- read.csv("neuron_only_MAST_DEGs_BR_ST.csv", row.names = 1)

head(degs3)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs3), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs3 <- degs3 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs3$ENTREZID, degs3$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=15, font.size=10) +
  ggtitle("GO Biological Process Enrichment Bioreactor vs. Static Neurons Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_BRvsST.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nBioreactor vs. Static Neurons Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_BRvsST.pdf", p, width = 8, height = 7)







##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################



# IPCs



#### Bioreactor vs MF IP


# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur
head(seur)
table(seur$sample)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712
unique(seur$predicted_CellType)

seur_BR_MF <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Microfluidic'))
seur_BR_MF_IP <- subset(seur_BR_MF, subset = predicted_CellType %in% c( 'IP'))
table(seur_BR_MF_IP$sample.simple)
table(seur_BR_MF_IP$predicted_CellType)

Idents(seur_BR_MF_IP) <- "sample.simple"

# Perform DEG analysis between conditions
degs3 <- FindMarkers(
  object = seur_BR_MF_IP,
  ident.1 = "Bioreactor",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs3 <- degs3 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs3 <- degs3 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs3, "IPC_only_MAST_DEGs_BR_MF.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs3 <- read.csv("neuron_only_MAST_DEGs_BR_ST.csv", row.names = 1)

head(degs3)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs3), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs3 <- degs3 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs3$ENTREZID, degs3$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=10, font.size=10) +
  ggtitle("GO Biological Process Enrichment Bioreactor vs. Microlfuidic IPC Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_BRvsMF_IPC.csv", row.names = FALSE)
getwd()

library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nBioreactor vs. Microfluidic IPC Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_BRvsMF_allIPC_downsampled.pdf", p, width = 8, height = 7)
getwd()
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################


#### Static vs MF IP


# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712
unique(seur$predicted_CellType)

seur_ST_MF <- subset(seur, subset = sample.simple %in% c( "Static", 'Microfluidic'))
seur_ST_MF_IP <- subset(seur_ST_MF, subset = predicted_CellType %in% c( 'IP'))
table(seur_ST_MF_IP$sample.simple)
table(seur_ST_MF_IP$predicted_CellType)

Idents(seur_ST_MF_IP) <- "sample.simple"

# Perform DEG analysis between conditions
degs3 <- FindMarkers(
  object = seur_ST_MF_IP,
  ident.1 = "Static",          
  ident.2 = "Microfluidic",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs3 <- degs3 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs3 <- degs3 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs3, "IPC_only_MAST_DEGs_ST_MF.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs3 <- read.csv("neuron_only_MAST_DEGs_BR_ST.csv", row.names = 1)

head(degs3)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs3), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs3 <- degs3 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs3$ENTREZID, degs3$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=10, font.size=10) +
  ggtitle("GO Biological Process Enrichment Static vs. Microlfuidic IPC Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_STvsMF_IPC.csv", row.names = FALSE)
getwd()

library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nStatic vs. Microfluidic IPC Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_STvsMF_allPC_downsampled.pdf", p, width = 15, height = 14)



##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

############ BR vs ST IPC

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712
unique(seur$predicted_CellType)

seur_BR_ST <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Static'))
seur_BR_ST_IP <- subset(seur_BR_ST, subset = predicted_CellType %in% c( 'IP'))
table(seur_BR_ST_IP$sample.simple)
table(seur_BR_ST_IP$predicted_CellType)

Idents(seur_BR_ST_IP) <- "sample.simple"

# Perform DEG analysis between conditions
degs3 <- FindMarkers(
  object = seur_BR_ST_IP,
  ident.1 = "Bioreactor",          
  ident.2 = "Static",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs3 <- degs3 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs3 <- degs3 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs3, "IPC_only_MAST_DEGs_BR_ST.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs3 <- read.csv("neuron_only_MAST_DEGs_BR_ST.csv", row.names = 1)

head(degs3)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs3), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs3 <- degs3 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs3$ENTREZID, degs3$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=10, font.size=10) +
  ggtitle("GO Biological Process Enrichment Static vs. Microlfuidic IPC Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_STvsMF_IPC.csv", row.names = FALSE)
getwd()

library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nBioreactor vs. Static IPC Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_BRvsST_allIP_downsampled.pdf", p, width = 8, height = 7)





##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

# BR vs. MF #need to update this block of code, started to work on the loop in the meantime 050425


# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(ggrepel)

setwd("~/Microfluidics McCain Collab/single_cell_run")
seur <- readRDS('Final_seur_quadrato_removed_unknown_downsampled_labelTransfer.Robj')

seur
head(seur)
# Define patterns for mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)

# Combine to a single vector of genes to remove
genes_to_remove <- c(mt_genes, rb_genes)

# Remove these genes from the Seurat object
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))

# Confirm removal
length(mt_genes)  # how many mitochondrial genes removed
length(rb_genes)  # how many ribosomal genes removed
length(genes_to_remove)  # total genes removed
dim(seur)  # new dimension of Seurat object
# [1] 44085 20712
unique(seur$predicted_CellType)

seur_BR_ST <- subset(seur, subset = sample.simple %in% c( "Bioreactor", 'Microfluidic'))
seur_BR_ST_IP <- subset(seur_BR_ST, subset = predicted_CellType %in% c( 'IP'))
table(seur_BR_ST_IP$sample.simple)
table(seur_BR_ST_IP$predicted_CellType)

Idents(seur_BR_ST_IP) <- "sample.simple"

# Perform DEG analysis between conditions
degs3 <- FindMarkers(
  object = seur_BR_ST_IP,
  ident.1 = "Bioreactor",          
  ident.2 = "Static",             
  test.use = "MAST",               # Use MAST for better handling of zero-inflated data
  latent.vars = c("nCount_RNA"),   # Control for library size or other covariates
  logfc.threshold = 0.25,          
  min.pct = 0.1                    
)


# Add a column to classify genes based on significance
degs3 <- degs3 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

degs3 <- degs3 %>%
  dplyr::mutate(significance = factor(significance, levels = c("Upregulated", "Downregulated"))) %>%
  dplyr::arrange(significance, p_val)

#output_dir <- "MAST_DEGs"
#dir.create(output_dir, showWarnings = FALSE)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
write.csv(degs3, "IPC_only_MAST_DEGs_BR_ST.csv", row.names = TRUE)

###################################################################################################################################
########## GO BP as in Paulsen et al., ##########################################

degs3 <- read.csv("neuron_only_MAST_DEGs_BR_ST.csv", row.names = 1)

head(degs3)

# Map gene symbols (rownames) to Entrez IDs
gene_df <- bitr(rownames(degs3), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Merge Entrez IDs into DEG dataframe
degs3 <- degs3 %>%
  mutate(gene = rownames(.)) %>%
  left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

# Split Entrez IDs by significance
gene_list <- split(degs3$ENTREZID, degs3$significance)

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
head(ego_compare)

# Plot results
dotplot(ego_compare, showCategory=10, font.size=10) +
  ggtitle("GO Biological Process Enrichment Static vs. Microlfuidic IPC Only")

# Save results
dir.create("output", showWarnings = FALSE)

write.csv(as.data.frame(ego_compare), "output/GO_BP_enrichment_results_STvsMF_IPC.csv", row.names = FALSE)
getwd()

library(ggplot2)
library(dplyr)
library(forcats)

# 1️⃣ Extract your compareCluster result
df <- as.data.frame(ego_compare)

# 2️⃣ Parse GeneRatio into numeric fraction
df <- df %>%
  mutate(
    GeneRatioNum = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                     function(z) as.numeric(z[1]) / as.numeric(z[2])))
  )

# 3️⃣ Top 10 per cluster
df_up   <- df %>% filter(Cluster=="Upregulated")   %>% arrange(p.adjust) %>% slice(1:10)
df_down <- df %>% filter(Cluster=="Downregulated") %>% arrange(p.adjust) %>% slice(1:10)

df_top <- bind_rows(df_up, df_down) %>%
  # assign custom x positions
  mutate(
    xpos        = if_else(Cluster=="Upregulated",   0.2, 0.4),
    Description = factor(Description, levels = rev(unique(Description)))
  )

# 4️⃣ Plot
p <- ggplot(df_top, aes(x = xpos, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  
  # custom x-axis
  scale_x_continuous(
    breaks = c(0.2, 0.4),
    labels = c("Upregulated", "Downregulated"),
    limits = c(0.1, 0.5)       # give a little padding
  ) +
  
  # black-to-light-grey for p.adjust
  scale_color_gradient(
    low  = "black",
    high = "grey80",
    name = "Adj. p-value"
  ) +
  
  # dot size = gene ratio
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Ratio"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment\nBioreactor vs. Static IPC Only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

# 5️⃣ Print & Save
print(p)
setwd("~/Microfluidics McCain Collab/single_cell_run/MAST_DEGs")
ggsave("results/GO_BP_Top10_UpDown_BRvsST_allIP_downsampled.pdf", p, width = 8, height = 7)







########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


# The loop!!!!


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
library(MAST)
library(Seurat)
library(forcats)
library(ggrepel)
library(glue)


setwd("J:/MF_Single_Cell_Updated_No_Static")
# Create output folder
output_dir <- "MAST_DEG_GO_top_5"


run_deg_go_loop <- function(seur,
                            celltypes,
                            output_dir = "MAST_DEG_GO_top_5",
                            conditions = c("Bioreactor", "Microfluidic"),
                            logfc_cutoff = 0.5,
                            min_cells = 50,
                            min_genes_required = 10) {
  
  # Load required packages
  require(dplyr)
  require(ggplot2)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(MAST)
  require(Seurat)
  require(glue)
  
  dir.create(output_dir, showWarnings = FALSE)
  
  for (ct in celltypes) {
    for (cmp in combn(conditions, 2, simplify = FALSE)) {
      cond1 <- cmp[1]
      cond2 <- cmp[2]
      contrast_label <- paste0(cond1, "_vs_", cond2)
      
      # Sanitize file-safe labels
      ct_safe <- gsub("[/\\: ]+", "_", ct)
      contrast_safe <- gsub("[/\\: ]+", "_", contrast_label)
      
      # Subset Seurat
      sub <- subset(seur, subset = sample.simple %in% c(cond1, cond2) & predicted_CellType == ct)
      if (ncol(sub) < min_cells) next
      
      Idents(sub) <- "sample.simple"
      message(glue::glue("Running DEG for {ct}: {cond1} vs {cond2}"))
      
      # DEG
      degs <- tryCatch({
        FindMarkers(
          object = sub,
          ident.1 = cond1,
          ident.2 = cond2,
          test.use = "MAST",
          latent.vars = c("nCount_RNA"),
          logfc.threshold = 0.25,
          min.pct = 0.1
        )
      }, error = function(e) {
        message(glue::glue("DEG failed for {ct} {contrast_label}: {e$message}"))
        return(NULL)
      })
      
      if (is.null(degs) || nrow(degs) == 0) next
      
      degs <- degs %>%
        mutate(gene = rownames(.)) %>%
        filter(p_val_adj < 0.05 & abs(avg_log2FC) >= logfc_cutoff) %>%
        mutate(
          significance = case_when(
            avg_log2FC >  logfc_cutoff ~ "Upregulated",
            avg_log2FC < -logfc_cutoff ~ "Downregulated"
          ),
          significance = factor(significance, levels = c("Upregulated", "Downregulated"))
        )
      
      if (nrow(degs) == 0 || length(unique(degs$significance)) == 0) next
      
      # Gene count filter before GO
      gene_counts <- degs %>%
        group_by(significance) %>%
        summarise(n = n(), .groups = "drop")
      too_few <- any(gene_counts$n[gene_counts$significance %in% c("Upregulated", "Downregulated")] < min_genes_required)
      if (too_few) {
        message(glue::glue("Skipping GO for {ct} {contrast_label}: too few genes in one or more groups"))
        next
      }
      
      # Save DEGs
      write.csv(degs, file = file.path(output_dir, paste0("DEG_", ct_safe, "_", contrast_safe, ".csv")), row.names = TRUE)
      
      # Entrez mapping
      gene_df <- bitr(degs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      degs <- left_join(degs, gene_df, by = c("gene" = "SYMBOL")) %>% filter(!is.na(ENTREZID))
      gene_list <- split(degs$ENTREZID, degs$significance)
      if (length(gene_list) == 0) next
      
      # GO:BP enrichment
      ego <- compareCluster(
        geneCluster = gene_list,
        fun = "enrichGO",
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      ego_df <- as.data.frame(ego)
      if (nrow(ego_df) == 0) next
      
      write.csv(ego_df, file = file.path(output_dir, paste0("GO_BP_", ct_safe, "_", contrast_safe, ".csv")), row.names = FALSE)
      
      # …after you compute `ego_df` and write out CSV…
      
      # Process GeneRatio
      ego_df <- ego_df %>%
        mutate(
          GeneRatioNum = as.numeric(sapply(
            strsplit(as.character(GeneRatio), "/"),
            function(z) as.numeric(z[1]) / as.numeric(z[2])
          ))
        )
      
      # Top GO terms (per direction)
      top_up   <- ego_df %>%
        filter(Cluster == "Upregulated") %>%
        arrange(p.adjust) %>%
        slice_head(n = 5)
      
      top_down <- ego_df %>%
        filter(Cluster == "Downregulated") %>%
        arrange(p.adjust) %>%
        slice_head(n = 5)
      
      df_top <- bind_rows(top_up, top_down) %>%
        mutate(
          # old xpos is no longer needed (you used 0.2/0.4 earlier)
          # Now create xpos2 purely from the Cluster label:
          xpos2 = case_when(
            Cluster == "Upregulated"   ~ 0.10,
            Cluster == "Downregulated" ~ 0.15
          ),
          # Preserve the ordering of the y‐axis (so the “top 5” show nicely):
          Description = factor(Description, levels = rev(unique(Description)))
        )
      
      # Dot plot
      p <- ggplot(df_top, aes(x = xpos2, y = Description)) +
        geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
        
        # keep the two categories at 0.10 and 0.20
        scale_x_continuous(
          breaks = c(0.10, 0.15),
          labels = c("Upregulated", "Downregulated")
        ) +
        
        # zoom in around those two x values—nothing is dropped
        coord_cartesian(xlim = c(0.09, 0.16)) +
        
        scale_color_gradient(
          low  = "black",
          high = "grey80",
          name = "Adj. p-value"
        ) +
        
        scale_size_continuous(
          range = c(3, 8),
          name = "Gene Ratio"
        ) +
        
        labs(
          title = glue::glue("GO BP Enrichment: {ct} — {contrast_label}"),
          x = NULL,
          y = NULL
        ) +
        
        theme_minimal(base_size = 14) +
        theme(
          plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          legend.title     = element_text(size = 12),
          legend.text      = element_text(size = 11),
          panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          
          # tighten margins and vertical spacing
          plot.margin      = margin(t = 5, r = 5, b = 5, l = 5),
          panel.spacing.y  = unit(0.1, "lines")
        )
      
      # Save the PDF
      ggsave(
        filename = file.path(
          output_dir,
          paste0("GO_BP_Top5_", ct_safe, "_", contrast_safe, ".pdf")
        ),
        plot   = p,
        width  = 8,
        height = 6,
        dpi    = 300
      )
      
    }
  }
}


########### RUN THE FUNCTION ########################
##############################



setwd("J:/MF_Single_Cell_Updated_No_Static")
seur <- readRDS('Final_seur_quadrato_MFandBioOnly_removedUnknown_downsampled_labelTransfer.Robj')

getwd()
# Remove mitochondrial and ribosomal genes
mt_genes <- grep("^MT-", rownames(seur), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(seur), value = TRUE)
genes_to_remove <- c(mt_genes, rb_genes)
seur <- subset(seur, features = setdiff(rownames(seur), genes_to_remove))


unique(seur@meta.data$predicted_CellType)

# Define conditions and cell types
celltypes_to_test <- c( 'oRG', 'IP', 'CFuPN', 
                        'CPN')  
output_dir <- "MAST_DEG_GO_top_5"

run_deg_go_loop(seur, celltypes = celltypes_to_test, output_dir = "MAST_DEG_GO", logfc_cutoff = 0.5,
                min_cells = 50,
                min_genes_required = 10
)
 




