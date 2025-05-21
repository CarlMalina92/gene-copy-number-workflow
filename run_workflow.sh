#!/bin/bash

# Exit on error
set -eu

# Define paths and variables
bam=~/Documents/HoppyYeast/HY123/ploidyAnalysis/test_data_gene_cn_workflow/HY123-sorted-rmdups.bam
reference=~/Documents/HoppyYeast/HY123/ploidyAnalysis/test_data_gene_cn_workflow/S288C_R64-3-1-chr_renamed.fasta
gff=~/Documents/HoppyYeast/HY123/ploidyAnalysis/test_data_gene_cn_workflow/S288C_R64-3-1-chr_renamed.gff
dict=~/Documents/HoppyYeast/HY123/ploidyAnalysis/test_data_gene_cn_workflow/S288C_R64-3-1-chr_renamed.dict
gatk=/Users/carlsberg/Documents/Software/gatk-4.6.2.0/gatk
outdir=~/Documents/HoppyYeast/HY123/ploidyAnalysis/test_data_gene_cn_workflow/results/
expected_ploidy=4.0
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