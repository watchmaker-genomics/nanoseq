#!/bin/bash

# Define the downsampling depths
DOWNSAMPLES=(2000 5000)

for N in "${DOWNSAMPLES[@]}"
do
    RUN_NAME="CM_test_downsample_${N}"
    nextflow run main.nf \
      --input /Users/campbell.mcduling/WMG_repos/research_analysis/GH1455/data/leo_ont_test_nanoseq_samplesheet.csv \
      --protocol cDNA \
      --skip_demultiplexing \
      --skip_fusion_analysis \
      --skip_differential_analysis \
      -profile awsbatch,docker -c aws_batch.config \
      --run ${RUN_NAME} \
      --downsample_depth ${N} \
      --skip_modification_analysis \
      --skip_xpore \
      --skip_m6anet \
      --quantification_method stringtie2
done
