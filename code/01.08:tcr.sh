#!/usr/bin/bash

# This script process the single-cell TCR data using the clustering output

REPORT=/home/ciro/ad_hoc/ucolitis/results/tcr/runs25to28_clean1_ccreg
SOBJECT=/mnt/beegfs/ciro/runs25to28_clean1_ccreg_object_lock_mean0.01_pct20_pc20.rds
AGGGEX=/home/ciro/ad_hoc/ucolitis/raw/cellranger/aggr/runs25to28/outs/aggregation.csv
TCRCODE=/home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/preliminary_TCR_data_analysis.1.5.4.R
TCRCODE=/home/ciro/scripts/crtcr/analysis.R
TAGS="c('origlib', 'celltype', 'orig.group', 'RNA_snn_res.0.6', 'orig.condition')"

# Build the input CSV for to aggregate
mkdir --parents ${REPORT}
AGGR=${REPORT}/aggregation.csv
echo "library_id,clonotypes,annotations,sample.sffx" > ${AGGR}
unset FNAMES; FNAMES=(`tail -n +2 ${AGGGEX} | cut -d, -f1 | sed -E 's/_Gex|^[0-9]{1,}_//g'`)
for I in `seq 0 $((${#FNAMES[@]}-1))`; do
  SNAME=${FNAMES[${I}]}; PNAME=$(ls -d /home/ciro/ad_hoc/ucolitis/raw/cellranger/vdj/*TCR | grep ${SNAME})
  echo ${SNAME}
  CLONOF="${PNAME}/outs/clonotypes.csv"; if [[ ! -f ${CLONOF} ]]; then echo "Missing ${CLONOF}"; fi
  ANNOTF="${PNAME}/outs/filtered_contig_annotations.csv"; if [[ ! -f ${ANNOTF} ]]; then echo "Missing ${ANNOTF}"; fi
  echo "${SNAME},${CLONOF},${ANNOTF},$((${I}+1))" >> ${AGGR}
done
cat ${AGGR}

# Aggregate clonotypes.csv and filtered_contig_annotations.csv
Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/aggr_vdj.2.2.R \
  --ReportsPath ${REPORT} \
  --GenInputPath ${REPORT} \
  --AggrTable ${AGGR}\
  --SampleCountThold 2 --FreqThold 2
ls -loh ${REPORT}

# Merge clonotypes into the Seurat metadata and create basic visualizations
Rscript ${TCRCODE} \
  --ReportsPath=${REPORT} \
  --TCRContigs=${REPORT}/filtered_contig_annotations_aggr.csv \
  --TCRClonotypes=${REPORT}/clonotypes_aggr.csv \
  --SeuratObj=${SOBJECT} \
  --Tags="${TAGS}"
ls -loh ${REPORT}
