#!/bin/bash

# Exit on error
set -eu

# Define paths and variables
bam=/path/to/bam
reference=/path/to/fasta_reference
gff=/path/to/gff
dict=/path/to/dict
gatk=/path/to/gatk_executable
outdir=/path/to/output_directory
expected_ploidy=
highlight_genes=

# Run the pipeline
gene-copy-number \
    --bam "$bam" \
    --reference "$reference" \
    --gff "$gff" \
    --dict "$dict" \
    --coverage-summary "${outdir}sample_gene_summary.txt" \
    --output-dir "$outdir" \
    --highlight-genes "$highlight_genes" \
    --create-dict \
    --gatk-path "$gatk" \
    --run-depth \
    --expected-ploidy "$expected_ploidy"