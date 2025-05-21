# Gene Copy Number Workflow

This repository contains a bioinformatics workflow for estimating gene copy number variations from aligned short-read sequencing data.

## Workflow Overview

This section outlines each step of the gene copy number estimation process.

### 1. Preprocessing

Prepares input data for analysis.

- **Inputs**: Aligned reads (BAM), reference genome (FASTA), gene annotations (GFF3).
- **Steps**:
  - Create a sequence dictionary.
  - Convert GFF3 annotations to an interval list.
  - Extract gene positions for analysis.

### 2. Coverage Analysis

Core computation for estimating copy numbers.

- Run GATK `DepthOfCoverage` or similar.
- Normalize read coverage.
- Estimate gene copy numbers.
- Determine statistical significance of deviations.

### 3. Visualization

Generates plots and summary outputs.

- Genome-wide copy number plots.
- Highlight genes of interest.
- Export significant gene lists.

### Outputs

- A `.tsv` file containing gene name, estimated copy number, z-score, p-value, adjusted p-value, and genomic position 
- A plot of genome-wide estimated copy numbers, highlighting significance
- Reproducible intermediate files

### Significance Detection (Statistical Overview)

1. **Copy Number Calculation**
   
  Copy number is estimated from normalized read coverage using the formula:
  ```
  estimated_cnᵢ = (mean_coverageᵢ / median_coverage) × expected_ploidy
  ```

  Where:

    - `mean_coverageᵢ` is the coverage for gene *i*
    - `median_coverage` is the median coverage across all genes
    - `expected_ploidy` is typically 2.0 for diploid genomes

2. **Outlier Filtering**
  
  Genes with extreme copy number are excluded before computing statistics, based on:
  ```
  IQR = Q3 - Q1
  Lower Bound = Q1 - 1.5 × IQR
  Upper Bound = Q3 + 1.5 × IQR
  ``` 

3. **Z-score calculation**
  
  Standard z-scores for gene copy number are computed using:
  ```
  zᵢ = (estimated_cnᵢ - μ) / σ
  ```

  Where:

    - μ is the mean estimated copy number across non-outlier genes
    - σ is the standard deviation of estimated copy numbers

4. **Two-sided p-value from z-score**
  
  Computes the probability of observing a deviation this extreme or more.

  ```
  pᵢ = 2 × (1 - Φ(|zᵢ|))
  ```
  Where:

    - Φ is the cumulative distribution function (CDF) of the standard normal distribution

5. **FDR Correction (Benjamini-Hochberg)**
  
  Adjusts for multiple testing:

  ```
  p_adj(i) = min((pᵢ × n) / rank(pᵢ), 1)
  ```

  Where:

    - `n` is the number of genes tested
    - `rank(pᵢ)` is the rank of pᵢ among all p-values (ascending)

## Installation

### Using Conda (Recommended)

1. Clone the repository:

```bash
git clone https://github.com/your-username/gene-copy-number-workflow.git
cd gene-copy-number-workflow
```

2. Create the conda environment:

```bash
conda env create -f environment.yml
```

3. Activate the environment:

```bash
conda activate gene-copy-number
```

4. Install the package locally:

```bash
pip install .
```

> **Note:** Even though `conda` installs all required dependencies, you still need to run:
>
> ```bash
> pip install .
> ```
>
> This installs your local workflow code and registers the `gene-copy-number` command-line tool. Without this step, the CLI won't be available.

## 🔧 Automatic File Indexing

The workflow automatically handles the creation of required reference and alignment indexes:

- **FASTA `.dict` file**: If missing and `--create-dict` is passed, created using:
  ```
  gatk CreateSequenceDictionary -R reference.fasta -O reference.dict
  ```

- **FASTA `.fai` file**: If missing, created using:
  ```
  samtools faidx reference.fasta
  ```

- **BAM `.bai` file**: If missing, created using:
  ```
  gatk BuildBamIndex -I sample.bam
  ```

## Testing the Workflow

After installation, you can test the workflow using the example data:

```bash
python tests/test_workflow.py
```

This test will:
- Use `test_data/sample_gene_summary.txt` and `test_data/gene_positions.tsv`
- Run `plot_genome_coverage()`
- Generate an output plot and check that it exists

You should see:
```
✅ Visualization test passed.
```

For manual testing, you can also inspect `test_data/test_output_plot.png`.

## Testing Data Preprocessing

You can validate the GFF and BAM processing logic using:

```bash
python tests/test_data_processing.py
```

This test will:

- Use `test_data/test_annotations.gff3` and `reference.dict`
- Require a BAM file at `test_data/test.bam` with a `.bai` index

If the BAM is missing, create it using:

```python
import pysam

header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}
with pysam.AlignmentFile("test_data/test.bam", "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "test_read"
    a.query_sequence = "ACTG" * 10
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 100
    a.mapping_quality = 60
    a.cigar = ((0, 40),)  # 40M
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0
    a.query_qualities = pysam.qualitystring_to_array("I" * 40)
    outf.write(a)

pysam.index("test_data/test.bam")
```
> **Note:** If you omit the `pysam.index()` step, the workflow will automatically use
> GATK to generate the BAM index by invoking:
>
> ```bash
> gatk BuildBamIndex -I test_data/test.bam
> ```
>
> This requires that `gatk` and Java ≥17 are correctly configured and passed via `--gatk-path`.

You should see:

```
✅ Data processing test passed.
```