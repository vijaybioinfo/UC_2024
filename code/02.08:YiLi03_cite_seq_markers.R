#!/usr/bin/R

####################
# CITE-seq markers #
####################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2022-10-20
# ---

# This script will create the table of markers from CITE-seq data

cat("**************************** CITE-seq data ****************************\n")

if(!exists("sc_adt_mat_sparse")){
  cat("Load in the ADT UMI matrix\n")
  fnames = paste0(list.files(
    path = mypath, pattern = "CITE", full.names = TRUE
  ), "/outs/raw_feature_bc_matrix")
  fnames = setNames(fnames, basename(gsub(".outs.*", "", fnames)))
  command = paste0("libnames = unique(", myobject, "@meta.data$origlib)")
  cat(command, "\n"); eval(parse(text = command))
  libnames = unique(sub("_Gex", "", libnames))
  fnames = fnames[grepl(paste0(libnames, collapse = "|"), fnames)]
  sc_adt_list = lapply(
    X = 1:length(fnames),
    FUN = function(x){
      cat(names(fnames)[x], "\n")
      y <- Seurat::Read10X(data.dir = fnames[[x]], strip.suffix = TRUE)
      colnames(y) <- paste0(colnames(y), "-", names(fnames)[x])
      y#[grep("MKR", rownames(y), value = TRUE), ]
    }
  ); sapply(sc_adt_list, function(x) head(rownames(x)) )
  i = 1; sc_adt_mat = sc_adt_list[[i]] # merging
  for(i in seq(i+1, length(sc_adt_list))){
    sc_adt_mat <- cbind(sc_adt_mat, sc_adt_list[[i]])
  }; # rm(sc_adt_list)
  class(sc_adt_mat)
  bc_maxs <- matrixStats::colMaxs(as.matrix(sc_adt_mat))
  bc_sums <- Matrix::colSums(sc_adt_mat); mean(bc_sums > 0); sum(bc_sums > 0)
  if(1){
    cat("MAX UMI -------\n% > 10:", mean(bc_maxs > 10), "\n");
    cat(" Total:", sum(bc_maxs > 10), "\n")
    cat("SUM UMI -------\n% > 10:", mean(bc_sums > 10), "\n");
    cat(" Total:", sum(bc_sums > 10), "\n")
    cat("\nSums > 0 & Max > 10", mean(bc_sums > 0 & bc_maxs > 10), "\n");
    cat(" Total:", sum(bc_sums > 0 & bc_maxs > 10), "\n")
  }
  sc_adt_mat <- sc_adt_mat[, bc_sums > 0 & bc_maxs > 10]
  rownames(sc_adt_mat) <- gsub("MKR\\-", "", rownames(sc_adt_mat))
  command = paste0("sco_cells = gsub('\\\\-.*', '-', rownames(", myobject, "@meta.data))")
  cat(command, "\n"); eval(parse(text = command))
  command = paste0("sco_cells = paste0(sco_cells,",
    "gsub('^[0-9]{1,}_|_Gex', '', ", myobject, "@meta.data$origlib))")
  cat(command, "\n"); eval(parse(text = command))
  colnames(sc_adt_mat) = gsub("(.*\\-)[0-9]{1,}_|_CITE", "\\1", colnames(sc_adt_mat))
  head(colnames(sc_adt_mat)); head(sco_cells) # matching libraries
  mean(colnames(sc_adt_mat) %in% sco_cells) # overlaping with expression data
  mean(sco_cells %in% colnames(sc_adt_mat)) # overlaping with CITE data
  command = paste0("tvar <- table(", myobject, "@meta.data$origlib, sco_cells %in% colnames(sc_adt_mat))")
  cat(command, "\n"); eval(parse(text = command))
  cbind(tvar, pct = tvar[, ncol(tvar)] / rowSums(tvar)); colSums(tvar)
  tvar <- sco_cells %in% colnames(sc_adt_mat)
  command = paste0("object_sub = ", myobject, "[, rownames(", myobject, "@meta.data)[tvar]]")
  cat(command, "\n"); eval(parse(text = command))
  sc_adt_mat_sparse = sc_adt_mat[, sco_cells[tvar]]; rm(sc_adt_mat)
  command = paste0("colnames(sc_adt_mat_sparse) = rownames(", myobject, "@meta.data)[tvar]")
  cat(command, "\n"); eval(parse(text = command))
  cat("Finished!\n")
  rm(sc_adt_list, fnames, bc_maxs, bc_sums, command)
}
