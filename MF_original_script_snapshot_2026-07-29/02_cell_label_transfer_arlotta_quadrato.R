
library(Seurat)

setwd("C:/Users/jeanp/Documents/Microfluidics McCain Collab/single_cell_run")

load('final_cca_object.Robj')
head(object.cca@meta.data)
object.cca

object.cca<- JoinLayers(object.cca)

object.cca<- NormalizeData(object.cca, normalization.method = 'LogNormalize', scale.factor = 10000)
object.cca<- FindVariableFeatures(object.cca)
object.cca <- ScaleData(object.cca)
object.cca <- RunPCA(object.cca, verbose = FALSE)
object.cca <- FindNeighbors(object.cca, dims = 1:30)
object.cca <- FindClusters(object.cca)
object.cca <- RunUMAP(object.cca, dims = 1:30)

# Ensure 'data.set' is a factor
object.cca$data.set <- as.factor(object.cca$data.set)

# Split into reference (Arlotta) and query (Quadrato)
ref_obj <- subset(object.cca, subset = data.set == "arlotta")
query_obj <- subset(object.cca, subset = data.set == "quadrato")


# 6) Transfer labels from reference
#    Use the same normalization.method and dims as above
transfer.anchors <- FindTransferAnchors(
  reference           = ref_obj,
  query               = query_obj,
  normalization.method = "LogNormalize",
  dims                  = 1:30
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata   = ref_obj$CellType,  # metadata column with labels
  dims      = 1:30
)

# 7) Add the predicted labels into the integrated object’s metadata; in our case the cca object
object.cca <- AddMetaData(object.cca, metadata = predictions)


head(object.cca)
object.cca

test.obj <- subset(object.cca, subset = data.set == "quadrato")

DimPlot(test.obj, group.by = 'predicted.id', cols = celltype_colors )

