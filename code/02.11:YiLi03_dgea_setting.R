#!/usr/bin/R

cat("====================== Setting DGEA comparisons ======================\n")

config_f = "/home/ciro/ucolitis/scripts/YiLi03_dgea_celltype-perturbation.yaml"
comp_f = sub("scripts", "info", sub("yaml", "csv", config_f))
config = yaml::read_yaml("/home/ciro/scripts/dgea/config.yaml")

if(1){
  config$method = "mastlog2cpm"
  config$project = "YiLi03_f1"
  config$metadata = "/home/ciro/ad_hoc_temp/ucolitis/results/clustering/YiLi03_f1/cite_seq/metadata.rds"
  config$expression_data = "/home/ciro/ad_hoc/ucolitis/raw/cellranger/aggr/mm3/outs"
  config$output_dir = "/home/ciro/ad_hoc_temp/ucolitis/results/dgea"
  config$comparisons <- config$comparisons[1] # only the file
  config$comparisons$file = NULL

  cat("---- WT vs KO per cell type\n") ## ------------------------------------------------
  temp = list("CD4", "CD8")
  for(i in 1:length(temp)){
    celltype = temp[[i]]
    cat(celltype, "\n")
    celltype_name = paste0(make.names(celltype), collapse = "-")
    context_i <- paste0("perturbation_", celltype_name)
    cat(context_i, "\n")
    config$comparisons[[context_i]] <- list(
      context = context_i, test_column = "orig.library_condition",
      contrast = c("WT", "KO"),
      filters = setNames(list(c(celltype)), "celltype_citeseq"),
      job = list(mem = "80gb")
    )
  }
  cat("Config:", config_f, "\n")
  yaml::write_yaml(config, file = config_f)
}
# Rscript ~/scripts/dgea/R/dgea_jobs.R -y /home/ciro/ucolitis/scripts/YiLi03_dgea_celltype-perturbation.yaml -s TRUE
