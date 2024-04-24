#!/usr/bin/bash

# ------------------------------------------------------------------------------
# title: YiLi03 mice data processing.
# purpose: This script was used to process our YiLi03 data from
#   Bcl6 KO vs WT CD4 T cells and CD8 T cells from colonic tissue from
#   mice with colitis.
#
# author:
#   - name: Ciro Ramírez-Suástegui
#     affiliation: La Jolla Institute for Immunology
#     email: ciro@lji.org, cramsuig@gmail.com
# Creation: 2022-08-02 Tue 13:55:53 PDT
# Last revised: 2022-10-19
#
# Copyright (C) 2022 Ciro Ramírez-Suástegui (cramsuig@gmail.com)
# Permission to copy and modify is granted under the GPL-3.0-or-later license
# ------------------------------------------------------------------------------

# You need:
# Sample/subject metadata
# Library metadata
# Feature barcodes sheet [optional]

PRJNAME=ucolitis
ADHOC=/mnt/BioAdHoc/Groups/vd-vijay/${USER}
ADHOCT=/mnt/bioadhoc-temp/Groups/vd-vijay/${USER}
DATASETID=YiLi03
MDATA_LIB=${HOME}/${PRJNAME}/info/YiLi_metadata_library.csv
FBARCODES=${HOME}/${PRJNAME}/info/fbarcodes_NV084.csv

PFILE=${HOME}/${PRJNAME}/scripts/${DATASETID}
EFILE=${ADHOC}/${PRJNAME}/raw/cellranger/aggr/mm3/outs/aggregation.csv
MDATA=${ADHOC}/${PRJNAME}/results/doublets/${DATASETID}/scrublet.rds

echo -e "\033[0;31m### Metadata ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\033[0m"
# # Locally
# scp ${HOME}/ciro\@lji.org\ -\ Google\ Drive/My\ Drive/ucolitis_ibd/metadata/{YiLi_metadata_library,fbarcodes_NV084}.csv \
#   ciro@herman-login2.liai.org:${HOME}/ucolitis/info
# In work machine
Rscript /home/ciro/scripts/functions/csvCorrect.R ${MDATA_LIB}
Rscript /home/ciro/scripts/functions/csvCorrect.R ${FBARCODES}

echo -e "\033[0;31m### Mapping ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\033[0m"
CRYAML0=https://raw.githubusercontent.com/vijaybioinfo/cellranger_wrappeR/e895b87916bce44018a52d2bbd736e53677fe732/config.yaml
CRYAMLT=$(mktemp).yaml

wget ${CRYAML0} -O ${CRYAMLT} --quiet
function join_by () {
  local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}";
}
cat ${CRYAMLT} |
  grep -vE "bcl2fastq|library_pattern| barcodes|sample_name1|additional" |
  sed 's| #.*||g;
  s|output_dir:.*|output_dir: '"${ADHOC}/${PRJNAME}"'/raw/cellranger|g;
  s|fastqs_dir:.*|fastqs_dir: '[\""$(join_by "\", \"" $(ls -d /mnt/BioAdHoc/Groups/vd-vijay/cramirez/seqteam/raw/{NV084,NV085}))"\"]'|g;
  s|fungal_allergy|'"${PRJNAME}"'|g;
  s|run:.*|run: '[\""$(join_by "\", \"" $(ls -d /mnt/NovaSeq/*{NV084,NV085}*))"\"]'|g;
  s|samples:.*|samples: '"${DATASETID}"'|g;
  s|hg19-3.0.0|mm10-3.0.0|g;
  s|main:.*|main: '"${FBARCODES}"'|g;
  s|GRCh38.*|GRCm38-alts-ensembl-3.1.0|g;
  s|submit:.*|submit: yes|g;
  s|aggregation:.*|aggregation: '"${MDATA_LIB}"'|g
' > ${PFILE}_cellranger.yaml
cat ${PFILE}_cellranger.yaml

sh /home/ciro/scripts/cellranger_wrappeR/run.sh -y ${PFILE}_cellranger.yaml -v

echo -e "\033[0;31m### Doublets/Scrublet ### %%%%%%%%%%%%%%%%%%%%%%%%%%%\033[0m"

conda activate doublets
cat /home/ciro/scripts/doublet_detection/config.yaml |
  grep -vE "class|dims" |
  sed 's|project_name:.*|project_name: '"${DATASETID}"'|;
  s|input_matrix:.*|input_matrix: '"${EFILE}"'|;
  s|metadata:.*|metadata: '"${MDATA_LIB}"'|;
  s|output_dir:.*|output_dir: '"${ADHOC}/${PRJNAME}"'/results/doublets|;
  s|score:.*|score: 0.4|;
' > ${PFILE}_scrublet.yaml
cat ${PFILE}_scrublet.yaml
/share/apps/R/3.6.1/bin/Rscript /home/ciro/scripts/doublet_detection/scrublet.R -y ${PFILE}_scrublet.yaml
conda deactivate

echo -e "\033[0;31m### Quality Control ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%\033[0m"
cat /home/ciro/scripts/quality_control/config.yaml |
  grep -vE "#.*|.*orig.HT_ID.global.*|.*subset.*" |
  sed 's|project_test|'"${DATASETID}"'_init|g;
  s|metadata:.*|metadata: '"${MDATA}"'|g; s|fungal_allergy|preethi|g;
  s|output_dir:.*|output_dir: '"${ADHOC}/${PRJNAME}"'/results/quality_control|;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g
' > ${PFILE}_qc0.yaml
/share/apps/R/3.6.1/bin/Rscript /home/ciro/scripts/quality_control/single_cell.R -y ${PFILE}_qc0.yaml

cat ${PFILE}_qc0.yaml |
  sed 's|project_name:.*|project_name: '"${DATASETID}"'_f1|g;
  s|nFeature_RNA:.*|nFeature_RNA: [500, 4000, 1]|g;
  s|nCount_RNA:.*|nCount_RNA: [100, 15000, 1]|g;
  s|percent.mt:.*|percent.mt: [-Inf, 4, 1]|g;
  /.*file.*/d;
' > ${PFILE}_qc1.yaml

if [[ "${DATASETID}" == "example" ]]; then
  sed -i '$ d' ${PFILE}_qc1.yaml # '$ d' remove last line
  echo -e "  expr: \"(orig.location == 'LNT' & nFeature_RNA >= 1500) | (orig.location == 'TT' & nFeature_RNA >= 1000)\"\n..." >> ${PFILE}_qc1.yaml
fi
/share/apps/R/3.6.1/bin/Rscript /home/ciro/scripts/quality_control/single_cell.R -y ${PFILE}_qc1.yaml

echo -e "\033[0;31m### Clustering ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\033[0m"
QCLINE=$(/share/apps/R/3.6.1/bin/R --slave -e "
  x=yaml::read_yaml('${PFILE}_qc1.yaml')[['filtering']]; xx=x[c('nFeature_RNA', 'nCount_RNA', 'percent.mt')];
  y=paste(sapply(names(xx), function(f) c(paste(f,'>=',x[[f]][1],collapse=' & '),paste(f,'<=',x[[f]][2],collapse=' & ')) ),collapse=' & ')
  if(isFALSE(is.null(x[['expr']]))) y=paste0('(', x[['expr']], ') & ', y)
  y=gsub('\\\&','\\\\\\\\&',y);y=gsub('\\\|','\\\\\\\\|',y);cat(y)
")
function lambda() { grep ${1} ${2} | sed 's/.*: //g'; }
MDATA_QC=$(lambda quality_control ${PFILE}_qc1.yaml)/$(lambda project_name ${PFILE}_qc1.yaml)/metadata_prefilter.rds

cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_QC}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|expr:.*|expr: "'"${QCLINE}"' \& doublet_scores <= '"$(lambda score ${PFILE}_scrublet.yaml)"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering0.yaml
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering0.yaml
# cd /mnt/bioadhoc-temp/Groups/vd-vijay/ciro/ucolitis/results/clustering/YiLi03_f1
# /share/apps/R/3.6.1/bin/Rscript /home/ciro/scripts/clustering/R/report.R --path ./ -c FALSE

cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1_cellcycle|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_QC}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|percent.mt]|percent.mt, cellcycle]|;
  s|expr:.*|expr: "'"${QCLINE}"' \& doublet_scores <= '"$(lambda score ${PFILE}_scrublet.yaml)"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering0_cellcycle.yaml
cat ${PFILE}_clustering0_cellcycle.yaml
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering0_cellcycle.yaml

cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1_wtkoreg|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_QC}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|percent.mt]|percent.mt, orig.library_condition]|;
  s|expr:.*|expr: "'"${QCLINE}"' \& doublet_scores <= '"$(lambda score ${PFILE}_scrublet.yaml)"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering0_wtkoreg.yaml
cat ${PFILE}_clustering0_wtkoreg.yaml
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering0_wtkoreg.yaml

CELLTYPE=CD8
CELLTYPE=CD4
MDATA_CITE=${ADHOCT}/${PRJNAME}/results/clustering/YiLi03_f1/cite_seq/metadata.rds
cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1_'"${CELLTYPE}"'|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_CITE}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|expr:.*|expr: "celltype_citeseq == '"'${CELLTYPE}'"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering1_${CELLTYPE}.yaml
cat ${PFILE}_clustering1_${CELLTYPE}.yaml # need to manually edit the filter
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering1_${CELLTYPE}.yaml

CELLTYPE=CD8
CELLTYPE=CD4
MDATA_CITE=${ADHOCT}/${PRJNAME}/results/clustering/YiLi03_f1/cite_seq/metadata.rds
cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1_'"${CELLTYPE}"'_wtkoreg|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_CITE}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|percent.mt]|percent.mt, orig.library_condition]|;
  s|expr:.*|expr: "celltype_citeseq == '"'${CELLTYPE}'"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering1_${CELLTYPE}.yaml
cat ${PFILE}_clustering1_${CELLTYPE}.yaml
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering1_${CELLTYPE}.yaml

cat /home/ciro/scripts/clustering/config.yaml |
  grep -vE ".*metadata_.*" |
  sed 's|clustering_test|'"${DATASETID}"'_f1_'"${CELLTYPE}"'_harmony|g;
  s|input_expression:.*|input_expression: '"${EFILE/aggregation*/}"'|g;
  s|metadata:.*|metadata: '"${MDATA_CITE}"'|g;
  s|output_dir:.*|output_dir: '"${ADHOCT}/${PRJNAME}"'/results/clustering/|;
  s|type: pca|type: harmony, batch: "orig.library_condition"|;
  s|expr:.*|expr: "celltype_citeseq == '"'${CELLTYPE}'"'"}|g;
  s|\./|/home/ciro/scripts/clustering/|g;
  s|exec:.*|exec: /share/apps/R/3.6.1/bin/Rscript|g;
  /^#.*/d
' > ${PFILE}_clustering1_${CELLTYPE}_harmony.yaml
cat ${PFILE}_clustering1_${CELLTYPE}_harmony.yaml
sh /home/ciro/scripts/clustering/run.sh -y ${PFILE}_clustering1_${CELLTYPE}_harmony.yaml
