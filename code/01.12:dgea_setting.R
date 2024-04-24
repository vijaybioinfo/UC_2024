#!/usr/bin/R

cat("====================== Setting DGEA comparisons ======================\n")

config_f = "/home/ciro/ucolitis/scripts/dgea_celltype-inflammation.yaml"
comp_f = sub("scripts", "info", sub("yaml", "csv", config_f))
config = yaml::read_yaml("/home/ciro/scripts/dgea/config.yaml")

if(1){
  config$method = "mastlog2cpm"
  config$project = "runs25to28_clean1_ccreg"
  config$metadata = "/home/ciro/ad_hoc/ucolitis/results/figures/metadata.rds"
  config$expression_data = "/home/ciro/large/ucolitis/results/clustering/runs25to28_clean1_ccreg/.object_stem_seurat_mean0.01_pct20.rds"
  config$output_dir = "/home/ciro/ad_hoc_temp/ucolitis/results/dgea"
  config$comparisons <- config$comparisons[1] # only the file
  config$comparisons$file = NULL

  cat("---- Inflammation\n") ## ------------------------------------------------
  temp = c(
    "0" = "CD4~T[FH]",
    "5" = "CD4~T[REF]",
    "7" = "CD4~T[REG]",
    "4" = "CD4~T[CM]^ANXA1",
    "6" = "CD4~T[CM]",
    "14" = "CD4~T[CM]^FOS",
    "10" = "CD8",
    "2" = "CD8~T^GZMK",
    "13" = "CD8~T[FH]",
    "1" = "CD8~T[RM]",
    "3" = "CD8~T[RM]^GZMA",
    "12" = "CD8~T[RM]^FCER1G",
    "9" = "CD8~T[RM]^XCL2",
    "8" = "CD8~T[H]*17",
    "11" = "CD4~T[H]*17",
    "15" = "Cycling"
  )
  temp = list("CD4~T[H]*17", "CD8~T[H]*17", "CD4~T[FH]", "CD8~T[FH]", c("CD4~T[REG]", "CD4~T[REF]"))
  for(i in 1:length(temp)){
    celltype = temp[[i]]
    cat(celltype, "\n")
    celltype_name = paste0(make.names(celltype), collapse = "-")
    context_i <- paste0("celltype_inflammation0_", celltype_name)
    cat(context_i, "\n")
    config$comparisons[[context_i]] <- list(
      context = context_i, test_column = "orig.group",
      contrast = c("Inflamed", "NonInvolved"),
      filters = setNames(list(c(celltype)), "celltype_subset"),
      job = list(mem = "40gb")
    )
    context_i <- paste0("celltype_inflammation1_", celltype_name)
    cat(context_i, "\n")
    config$comparisons[[context_i]] <- list(
      context = context_i, test_column = "orig.group",
      contrast = c("Inflamed", "NonInvolved~Uninflamed"),
      filters = setNames(list(c(celltype)), "celltype_subset"),
      job = list(mem = "40gb")
    )
  }
  cat("Config:", config_f, "\n")
  yaml::write_yaml(config, file = config_f)
}
# Rscript ~/scripts/dgea/R/dgea_jobs.R -y /home/ciro/ucolitis/scripts/dgea_celltype-inflammation.yaml -s TRUE
