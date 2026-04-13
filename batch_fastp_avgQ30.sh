#!/bin/bash

# Number of threads to use
THREADS=4
OUTDIR="fastp_filtered"
SUMMARY="${OUTDIR}/fastp_read_summary.csv"

# Create output directory
mkdir -p "$OUTDIR"

# Initialize summary file
echo "sample,total_reads_before,total_reads_after,reads_lost,percent_lost" > "$SUMMARY"

# Loop over all R1 FASTQ files
for r1 in *_R1.fq.gz; do
    base=$(basename "$r1" _R1.fq.gz)
    r2="${base}_R2.fq.gz"

    # Output filenames
    out_r1="${OUTDIR}/${base}_R1.filtered.fq.gz"
    out_r2="${OUTDIR}/${base}_R2.filtered.fq.gz"
    json="${OUTDIR}/${base}.json"
    html="${OUTDIR}/${base}.html"

    echo "Processing $base..."

    # Run fastp with average quality filtering
    fastp -i "$r1" -I "$r2" \
          -o "$out_r1" -O "$out_r2" \
          --average_qual 30 \
          --disable_adapter_trimming \
          --thread $THREADS \
          --json "$json" \
          --html "$html"

    # Parse counts from JSON report using jq
    before=$(jq '.summary.before_filtering.total_reads' "$json")
    after=$(jq '.summary.after_filtering.total_reads' "$json")
    lost=$((before - after))
    perc_lost=$(awk -v b="$before" -v a="$after" 'BEGIN { printf "%.2f", (b-a)/b*100 }')

    # Save to summary CSV
    echo "$base,$before,$after,$lost,$perc_lost" >> "$SUMMARY"

    echo "Done with $base"
    echo "----------------------------------"
done

echo "Summary saved to $SUMMARY"
