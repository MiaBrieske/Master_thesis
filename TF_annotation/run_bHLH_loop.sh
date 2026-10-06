#!/bin/bash


# Usage: ./run_bHLH_loop.sh 3
#        ./run_bHLH_loop.sh 4

BATCH="$1"

if [ -z "$BATCH" ]; then
    echo "Usage: $0 <batch_number>"
    exit 1
fi

cd /vol/analysis/bHLH\_annotator || exit 1

INPUT_DIR="/vol/analysis/plants/pep_batch_${BATCH}"
OUT_DIR="/vol/analysis/plants/bHLH_out"
LOG_DIR="${OUT_DIR}/log"

mkdir -p "$OUT_DIR" "$LOG_DIR"

for f in "${INPUT_DIR}"/*.fasta; do
	base=$(basename "$f" .fasta)

	python3 bHLH_annotator.py \
		--subject "$f" \
		--out "${OUT_DIR}/${base}_bHLH" \
		--info data/bHLH_annotator.csv \
		--mode_aln mafft \
		> "${LOG_DIR}/${base}_bHLH.log" 2>&1

done


