#!/usr/bin/R

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; redb = crayon::red$bold; green = crayon::green
}else{ cyan = redb = green = c }

{ cat(redb("### Locking objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/ciro/scripts/clustering/R/utilities.R")
  clustering_dir = "/home/ciro/ad_hoc_temp/ucolitis/results/clustering/"
  configs = list(
    list(name = "YiLi03_f1", sufix = "_seurat_mean0.01_pct25", pc = "15")
  )
  for (config in configs) {
    setwd(paste0(clustering_dir, config$name))
    mycells_f = paste0(
      "/mnt/beegfs/", Sys.info()[["user"]], "/", basename(getwd()), "_",
      gsub(".*seurat", "object_lock", config$sufix), "_pc", config$pc, ".rds")
    cat("File:", mycells_f, "\n")
    if(file.exists(mycells_f)) next
    mycells = complex_object_fetch(
      object_or_dir = paste0(".object_stem", config$sufix, ".rds"),
      id = paste0(config$sufix, "_pc", config$pc, ".rds"), verbose = TRUE
    )
    tvar <- grep("dim_", colnames(mycells@meta.data), inver = TRUE)
    mycells@meta.data <- mycells@meta.data[, tvar]
    names(mycells@reductions) <- gsub("_.*", "", names(mycells@reductions))

    cat("Saving object\n"); saveRDS(mycells, file = mycells_f);
    rm(mycells_f, mycells, config, tvar); gc()
  }; rm(configs)
}
