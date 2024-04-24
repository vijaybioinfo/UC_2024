#!/usr/bin/R

###############
# Global data #
###############

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2022-10-20
# ---

# This script loads all the files used in a global way

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; redb = crayon::red$bold; green = crayon::green
}else{ cyan = redb = green = c }

{ cat(redb("------------------ Setting global parameters -----------------\n"))
  cat(cyan("------------------ File names and initial objects\n"))
  fig_dir = "/home/ciro/ad_hoc_temp/ucolitis/results/clustering/YiLi03_f1"
  redu = list(umap = c("UMAP_1", "UMAP_2"))
  sc_tcells_clust = "RNA_snn_res.0.1"
  global_objects_f = c(
    gcolours = "/home/ciro/ucolitis/info/global_colours.csv",
    sc_tcells = "/mnt/beegfs/ciro/YiLi03_f1_object_lock_mean0.01_pct25_pc15.rds"
  )
  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R",
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/clustering/R/utilities.R",
    "/home/ciro/scripts/figease/figease.R",
    "cowplot", "patchwork", "tidyverse"
  )
  celltype_subset = c(
    "0" = "CD8~Tcell",
    "1" = "CD4~T[H]",
    "2" = "CD4~T[H]",
    "3" = "CD8~GZMA",
    "4" = "Cycling",
    "5" = "Macrophage",
    "6" = "CD4~T[REG]"
  )
  sc_tcells_ident = list(
    celltype = gsub("~.*", "", celltype_subset),
    celltype_subset = celltype_subset,
    order = names(celltype_subset)
  )
}
{ # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  if(!exists("include")) include = names(global_objects_f)
  dir.create(fig_dir); setwd(fig_dir)
  cat("Working at (fig_dir):", fig_dir, "\n")
  cat(cyan("------------------ Packages and functions\n")) # -------------------
  loaded <- lapply(
    X = packages_funcs,
    FUN = function(x){
      cat("*", x, "\n")
      if(!file.exists(x)){
        suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
      }else{ source(x) }
  }); theme_set(theme_cowplot())
  cat(cyan("------------------ Objects\n")) # ----------------------------------
  object_names = fig_global_objects(
    global_objects_files = global_objects_f,
    object_names = ls(pattern = "_f$|^redu|_clust$|_ident$"),
    reading_fun = 'readfile', stringsAsFactors = FALSE, check.name = FALSE, row.names = 1
  )
  sc_tcells_ident$colours = v2cols(names(sc_tcells_ident[[2]]), c(sc_tcells_ident$colours, gcolours))
  cat(green("------------------ Identities\n")) # ------------------------------
  ident_set_object(grep("sc_[[:alnum:]]{1,}$", names(global_objects_f), value = TRUE))
  print(format(
    object_names,
    justify = "centre", quote = FALSE)
  ); rm(include)
}
