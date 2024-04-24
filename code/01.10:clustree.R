#!/usr/bin/R

####################
# Clustering trees #
####################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-07-01
# ---

# This script allowed us to see what the cells were calssified as in each
# clustering analysis
# Used to detect cluster-name changes across these two clustering analyses

setwd("/home/ciro/large/ucolitis/results/clustering")
source("/home/ciro/scripts/clustering/R/clustree.R")
clustree_plot(
  dpath = c(
    "archive/2021-06-28/runs25to28_clean1_ccreg/.object_meta.data_seurat_mean0.01_pct20_pc20.rds",
    "runs25to28_clean1_ccreg_tmp/.object_meta.data_seurat_mean0.01_pct20_pc20.rds"
  ),
  column_cluster = c("snn_res.", "res.1$|0.9"), make_unique_columns = TRUE,
  sobject = "archive/2021-06-28/runs25to28_clean1_ccreg/object_lock_seurat_mean0.01_pct20_pc20.rds",
  output_dir = "runs25to28_clean1_ccreg/clustree_compare/"
)
