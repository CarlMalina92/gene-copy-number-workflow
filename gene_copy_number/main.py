import argparse
import os
import re
from gene_copy_number.data_processing import DataPreprocessor
from gene_copy_number.visualization import CopyNumberVisualizer

def main():
    parser = argparse.ArgumentParser(description="Gene Copy Number Estimation Workflow")

    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--gff", required=True, help="Gene annotations in GFF3 format")
    parser.add_argument("--dict", required=True, help="Sequence dictionary file (.dict)")
    parser.add_argument("--coverage-summary", required=True, help="Path to GATK DepthOfCoverage _gene_summary file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--highlight-genes", help="Comma-separated list of genes to highlight")
    parser.add_argument("--expected-ploidy", type=float, default=2.0, help="Expected ploidy (default: 2.0)")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR-adjusted p-value cut-off (default: 0.05)")
    parser.add_argument("--create-dict", action="store_true", help="Create sequence dictionary if not present")
    parser.add_argument("--gatk-path", default="gatk", help="Path to GATK executable (default: gatk)")
    parser.add_argument("--run-depth", action="store_true", help="Run GATK DepthOfCoverage automatically")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Check Java version
    try:
        import subprocess
        result = subprocess.run(["java", "-version"], capture_output=True, text=True)
        version_output = result.stderr.splitlines()[0] if result.stderr else result.stdout.splitlines()[0]
        match = re.search(r'version "(\d+)', version_output)
        if not match or int(match.group(1)) < 17:
            raise RuntimeError("Java 17 or later is required to run GATK. Please update your Java installation.")
        print(f"Detected Java version: {match.group(1)}")
    except Exception as e:
        raise RuntimeError(f"ava is not correctly installed or accessible: {e}")
    
    # Step 1: Initialize data processor and validate BAM file
    dp = DataPreprocessor(reference_genome=args.reference)
    dp.validate_bam_file(args.bam, gatk_path=args.gatk_path)

    # Step 2: Create a sequence dictionary (optional)
    if args.create_dict or not os.path.exists(args.dict):
        print("Creating sequence dictionary...")
        dp.create_sequence_dictionary(output_dict=args.dict, picard_path=args.gatk_path)

    # Step 3: Convert GFF3 to interval list
    interval_list = os.path.join(args.output_dir, "genes.interval_list")
    gene_positions = os.path.join(args.output_dir, "gene_positions.tsv")
    plot_output = os.path.join(args.output_dir, "copy_number_plot.png")

    count = dp.convert_gff_to_interval_list(
        gff_file = args.gff,
        dict_file = args.dict,
        output_file = interval_list,
        gene_positions_file=gene_positions
    )
    print(f"Converted {count} genes to interval list")

    # Step 4: Run GATK DepthOfCovergae (optional)
    if args.run_depth:
        depth_output_prefix = os.path.join(args.output_dir, "coverage")
        dp.run_depth_of_coverage(
            bam_path=args.bam,
            reference_path=args.reference,
            interval_list=interval_list,
            output_prefix=depth_output_prefix,
            gatk_path=args.gatk_path
        )
        args.coverage_summary = f"{depth_output_prefix}.sample_interval_summary"
    
    # Step 5: Plot copy number variation (CNV)
    if not args.coverage_summary or not os.path.exists(args.coverage_summary):
        raise FileNotFoundError("Coverage summary file is required and was not found.")
    
    visualizer = CopyNumberVisualizer()
    highlight = [gene.strip() for gene in args.highlight_genes.split(",")] if args.highlight_genes else None

    visualizer.plot_genome_coverage(
        coverage_file = args.coverage_summary,
        gene_positions_file = gene_positions,
        output_file = plot_output,
        expected_ploidy = args.expected_ploidy,
        highlight_genes = highlight,
        alpha = args.alpha
    )

if __name__ == "__main__":
    main()