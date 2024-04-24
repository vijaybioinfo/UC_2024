#!/usr/bin/R

# ------------------------------------------------------------------------------
# title: R Script Template.
# purpose: This is a template of best practices for an R script.
#
# author:
#   - name: Ciro Ramírez-Suástegui
#     affiliation: La Jolla Institute for Immunology
#     email: ciro@lji.org, cramsuig@gmail.com
# Creation: 2022-10-20 Thu 23:13:24 CEST
# Last revised: 2022-10-20
#
# Copyright (C) 2022 Ciro Ramírez-Suástegui (cramsuig@gmail.com)
# Permission to copy and modify is granted under the GPL-3.0-or-later license

# ------------------------------------------------------------------------------

source("/home/ciro/ucolitis/scripts/YiLi03_clustering0_lock.R")
source("/home/ciro/ucolitis/scripts/YiLi03_clustering0_global.R")
suppressPackageStartupMessages({
  library(Seurat)
})

{ cat(redb("### CITE-seq data ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  result_id = "cite_seq/"; dir.create(result_id, showWarnings = FALSE)
  mypath = "/home/ciro/ad_hoc/ucolitis/raw/cellranger/count"; myobject = "sc_tcells"
  source("/home/ciro/ucolitis/scripts/YiLi03_cite_seq_markers.R")
  object_sub[["ADT"]] <- Seurat::CreateAssayObject(counts = sc_adt_mat_sparse)
  Assays(object_sub); rownames(object_sub[["ADT"]])

  pp_cite = list()
  fconfigs = list(
    list(
      result_id = "cite_seq/log1p_", cite_pos = 3.5, features = c("CD4-C0001", "CD8B-C0230"),
      vars_scatter = c("celltype", "cluster_t"),
      adt_rna = setNames(c("Cd4", "Cd8b1"), c("CD4-C0001", "CD8B-C0230")),
      plot = TRUE
    ),
    list(
      result_id = "cite_seq/clr_", cite_pos = 0.5, features = c("CD4-C0001", "CD8B-C0230"),
      vars_scatter = c("celltype", "cluster_t"),
      adt_rna = setNames(c("Cd4", "Cd8b1"), c("CD4-C0001", "CD8B-C0230")),
      plot = TRUE
    )
  )

  object_sub$cluster_t = gsub("_", " ", object_sub$celltype_subset)
  for (fconfig in fconfigs[2]) {
    # Normalize ADT data
    DefaultAssay(object_sub) <- "ADT"
    if(isTRUE(grepl("log1p", fconfig$result_id))){
      object_sub@assays$ADT@data <- log1p(object_sub@assays$ADT@counts)
    }else{
      object_sub <- NormalizeData(
        object_sub, normalization.method = "CLR",
        margin = 1 # features (1) or cells (2)
      )
    }
    DefaultAssay(object_sub) <- "RNA"

    fname = paste0(fconfig$result_id, "adt_avg_heatmap")

    # Mean per cluster
    if(!"dplyr" %in% (.packages())) suppressPackageStartupMessages(library(dplyr))
    prots = rownames(object_sub@assays$ADT@data)
    adt_plot <- cbind(t(as.matrix(object_sub@assays$ADT@data)), object_sub@meta.data) %>%
      group_by(cluster_t) %>%
      summarize_at(.vars = prots, .funs = mean) %>% as.data.frame
    rownames(adt_plot) = NULL
    adt_plot = tibble::column_to_rownames(adt_plot, "cluster_t")
    ddf = reshape2::melt(cbind(Group = rownames(adt_plot), adt_plot))
    pp_cite[[fname]] = ggplot(ddf, aes(x = Group, y = variable)) +
      geom_tile(aes(fill = value)) +
      scale_fill_gradientn(colours = viridis::viridis(25)) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 0.8),
        axis.line = element_blank(), axis.ticks = element_blank()
      ) + scale_x_discrete(labels = function(l) parse(text=l)) +
      labs(fill = NULL, x = NULL, y = NULL)

    if(!isFALSE(fconfig$plot)) pdf(paste0(fname, ".pdf"), width = 10)
    # pp_cite[[fname]] = pheatmap::pheatmap(t(adt_plot), color = viridis::viridis(25), fontsize_row = 8, border_color = NA)
    if(!isFALSE(fconfig$plot)) print(pp_cite[[fname]]); graphics.off()

    for(i in fconfig$vars_scatter){
      fname = paste0(fconfig$result_id, "adt_scatter_", i); cat(fname, "\n")
      p <- FeatureScatter(object_sub,
        feature1 = paste0("adt_", fconfig$features[1]),
        feature2 = paste0("adt_", fconfig$features[2]), group.by = i) +
        scale_color_discrete(labels = function(l) parse(text=l)) +
        guides(colour = guide_legend(override.aes = list(size = 6), title = NULL))
      p <- plot_add_quadrants(p, list(fconfig$cite_pos, fconfig$cite_pos), type = "percent") + labs(title = NULL)
      tvar <- make_grid(table(object_sub[[i]])); tmp <- c(14,14)#plot_size(tvar)
      if(!isFALSE(fconfig$plot)){
        pdf(paste0(fname, ".pdf"), width = tmp[1], height = tmp[2])
        print(p); print(plot_add_densities(p, "colors"))
        print(p + facet_wrap(~colors, ncol = tvar[2]))
        graphics.off()
      }; pp_cite[[fname]] = p
    }

    # CITE/Expression comparisons with UMAP
    for(i in 1:length(fconfig$adt_rna)){
      g = names(fconfig$adt_rna)[i]; g_rna = fconfig$adt_rna[[i]]
      fname = paste0(fconfig$result_id, "umap_", g); cat(fname, "\n")
      p1 <- FeaturePlot(object_sub, g, cols = c("lightgrey", "#00a9bf")) +
        ggtitle(paste(g, "protein")) + labs(x = "Dim 1", y = "Dim 2")
      p2 <- FeaturePlot(object_sub, g_rna, cols = c("lightgrey", "#bf0000")) +
        ggtitle(paste(g_rna, "RNA")) + labs(x = "Dim 1", y = "Dim 2")
      if(!isFALSE(fconfig$plot)){
        pdf(paste0(fname, ".pdf"), width = 12)
        print(p1 | p2); graphics.off()
      }; pp_cite[[paste0(fname, "_prot")]] = p1; pp_cite[[paste0(fname, "_rna")]] = p2
    }

    p1 <- VlnPlot(object = object_sub, features = paste0("adt_", fconfig$features[1]),
      pt.size = 0.1, group.by = "cluster_t") +
      theme(legend.position = "none", axis.text.x = element_blank()) +
      labs(x = "Antibody Level")
    p2 <- VlnPlot(object = object_sub, features =  paste0("adt_", fconfig$features[2]),
      pt.size = 0.1, group.by = "cluster_t") +
      theme(legend.position = "none") +
      scale_x_discrete(labels = function(l) parse(text=l)) + labs(x = "Antibody Level")
    pdf(paste0(fconfig$result_id, "adt_violin.pdf"), width = 10)
    # print(plot_rm_layer(p1, "oint") / plot_rm_layer(p2, "oint"))
    print(p1 / p2)
    graphics.off()

    p1 <- FeaturePlot(object_sub, features = paste0("adt_", fconfig$features), cols = c("lightgrey", "#00a9bf"))
    pdf(paste0(fconfig$result_id, "adt_umap.pdf"), width = 12)
    print(p1)
    graphics.off()

    # Checking for double-positive
    cite_df = FetchData(object = object_sub, vars = c("cluster_t", paste0("adt_", fconfig$features)))
    cite_df[, -1] = sapply(
      X = names(cite_df[, -1]),
      FUN = function(x){
        if(is.numeric(cite_df[, x]))
          ifelse(cite_df[, x] > fconfig$cite_pos, gsub("adt_", "", x), "0") else cite_df[, x]
    }); tmp = reshape2::melt(table(cite_df)); colnames(tmp) <- make.names(colnames(tmp))
    tvar <- paste0(paste0("adt_", make.names(fconfig$features)), collapse = "+")
    cite_num = reshape2::dcast(
      data = tmp,
      formula = paste("cluster_t ~", tvar)
    )

    rowSums(cite_num[, -1])
    ddf_pct = as.matrix(cite_num[, -1]); rownames(ddf_pct) <- as.character(cite_num[, 1])
    ddf_pct = round(as.matrix(as.data.frame.matrix(prop.table(
      t(ddf_pct),
      margin = 2
    ))) * 100, 2)
    rowSums(ddf_pct); colSums(ddf_pct)

    palettebreaks = seq(from = min(ddf_pct), to = max(ddf_pct), by = 5)
    mypalette = colorRampPalette(colors = c("white", "red"), space = 'Lab')
    source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
    pdf(paste0(fconfig$result_id, "adt_heatmap.pdf"), width = 12, height = 7, onefile = FALSE)
    x = pheatmap(mat = ddf_pct, labels_col = sapply(colnames(ddf_pct), function(l) parse(text=l) ),
      angle_col = "45", scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
      color = mypalette(length(palettebreaks) - 1)
    )
    graphics.off()
  }
  object_sub$cluster_t = NULL
}

{ cat(redb("### Separating based on CITE-seq ### %%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  df2devide = FetchData(object_sub, vars = c("celltype", "celltype_subset", "CD4-C0001", "CD8B-C0230"))
  head(df2devide)
  df2devide$celltype_citeseq = NA
  temp = (df2devide$celltype == "CD4" | df2devide[["adt_CD4-C0001"]] > 0.5) & df2devide[["adt_CD8B-C0230"]] < 0.5
  df2devide$celltype_citeseq[temp] = "CD4"
  temp = (df2devide$celltype == "CD8" | df2devide[["adt_CD8B-C0230"]] > 0.5) & df2devide[["adt_CD4-C0001"]] < 0.5
  df2devide$celltype_citeseq[temp] = "CD8"
  temp = df2devide[["adt_CD8B-C0230"]] < 0.5 & df2devide[["adt_CD4-C0001"]] < 0.5
  df2devide$celltype_citeseq[temp] = NA
  temp = df2devide$celltype == "Macrophage"
  df2devide$celltype_citeseq[temp] = "Macrophage"
  table(df2devide$celltype_subset, df2devide$celltype_citeseq, useNA="always")
  object_sub@meta.data$celltype_citeseq = df2devide[rownames(object_sub@meta.data), ]$celltype_citeseq
  for(i in c("celltype_citeseq")){
    fname = paste0(fconfig$result_id, "adt_scatter_", i); cat(fname, "\n")
    p <- FeatureScatter(object_sub,
      feature1 = paste0("adt_", fconfig$features[1]),
      feature2 = paste0("adt_", fconfig$features[2]), group.by = i) +
      scale_color_discrete(labels = function(l) parse(text=l)) +
      guides(colour = guide_legend(override.aes = list(size = 6), title = NULL))
    p <- plot_add_quadrants(p, list(fconfig$cite_pos, fconfig$cite_pos), type = "percent") + labs(title = NULL)
    tvar <- make_grid(table(object_sub[[i]])); tmp <- c(14,14)#plot_size(tvar)
    if(!isFALSE(fconfig$plot)){
      pdf(paste0(fname, ".pdf"), width = tmp[1], height = tmp[2])
      print(p); print(plot_add_densities(p, "colors"))
      print(p + facet_wrap(~colors, ncol = tvar[2]))
      graphics.off()
    }; pp_cite[[fname]] = p
  }
  fname = paste0(fconfig$result_id, "adt_umap_", i); cat(fname, "\n")
  p <- Seurat::DimPlot(object_sub, group.by = i)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(p)
  graphics.off()
}

mdata = FetchData(
  object_sub,
  vars = grep("cluster|ADT|snn_res|celltype$|celltype_s", colnames(object_sub@meta.data), value = TRUE, invert = TRUE)
)
str(mdata)
saveRDS(mdata, file = "cite_seq/metadata.rds")
