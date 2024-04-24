#!/usr/bin/R

# This scripts creates the object used for downstream analyses

source("/home/ciro/scripts/clustering/R/utilities.R")
setwd("/home/ciro/large/ucolitis/results/clustering/runs25to28_clean1_ccreg")
sufix = "_seurat_mean0.01_pct20"; pc = 20
mycells = complex_object_fetch(
  object_or_dir = paste0(".object_stem", sufix, ".rds"),
  id = paste0(sufix, "_pc", pc, ".rds"), verbose = TRUE
)
tvar <- grep("dim_", colnames(mycells@meta.data), inver = TRUE)
mycells@meta.data <- mycells@meta.data[, tvar]
names(mycells@reductions) <- gsub("_.*", "", names(mycells@reductions))

# Cell types
cd4 = c("0", "5", "7", "4", "6", "14", "11")
cd8 = c("10", "2", "13", "1", "3", "12", "9", "8")
mycells = ident_set_names(
  object = mycells,
  ident = list(celltype = c(
    setNames(rep("CD4", length(cd4)), cd4),
    setNames(rep("CD8", length(cd8)), cd8), "15" = "Cycling"
  )), cluster_column = "RNA_snn_res.1"
)
table(mycells@meta.data[, c("cluster", "celltype")])

mycells_f = paste0(
  "/mnt/beegfs/", Sys.info()[["user"]], "/", basename(getwd()), "_",
  gsub(".*seurat", "object_lock", sufix), "_pc", pc)
cat("Saving to:", mycells_f, "\n")
saveRDS(mycells, file = paste0(mycells_f, ".rds"))

# # Getting ready to visualise it
# library(reticulate)
# anndata <- import("anndata", convert = FALSE)
# Seurat::DefaultAssay(mycells)
# edata = t(as.matrix(Seurat::GetAssayData(object = mycells)))
# adata <- anndata$AnnData(
#     X = edata, var = mycells@assays$RNA@meta.features,
#     obs = data.frame(mycells@meta.data, check.names = FALSE, stringsAsFactors = FALSE),
#     obsm  = list(
#       "X_umap" = Seurat::Embeddings(mycells, reduction = "umap"),
#       "X_pca" = Seurat::Embeddings(mycells, reduction = "pca")
#     )
# )
# afile = paste0(mycells_f, ".h5ad")
# if(file.exists(afile)) file.remove(afile)
# anndata$AnnData$write(adata, afile)
