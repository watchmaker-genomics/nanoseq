#!/bin/bash

# defaults
INPUT="./TSRT004_samplesheet_nanoseq.csv"
RUNNAME_PREFIX="TSRT004_downsample"
DOWNSAMPLES=(1000 5000 10000 20000 50000)

usage() {
  echo "Usage: $0 [-i input_csv] [-r runname_prefix] [-d comma_separated_depths]" 1>&2
  exit 1
}

while getopts ":i:r:d:h" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    r) RUNNAME_PREFIX="$OPTARG" ;;
    d)
      IFS=',' read -r -a DOWNSAMPLES <<< "$OPTARG"
      ;;
    h) usage ;;
    :) usage ;;
    \?) usage ;;
  esac
done
shift $((OPTIND - 1))

for N in "${DOWNSAMPLES[@]}"; do
  RUN_NAME="${RUNNAME_PREFIX}_${N}"
  nextflow run main.nf \
    --input "$INPUT" \
    --protocol cDNA \
    --skip_demultiplexing \
    --skip_fusion_analysis \
    -profile awsbatch,docker -c aws_batch.config \
    --run "${RUN_NAME}" \
    --downsample_depth "${N}" \
    --skip_modification_analysis \
    --skip_xpore \
    --skip_m6anet \
    --skip_differential_analysis \
    --quantification_method bambu
done
