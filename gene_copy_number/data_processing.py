import os
import subprocess
import gffutils
import pysam

class DataPreprocessor:
    """Handle data preprocessing for copy number analysis"""

    def __init__(self, reference_genome):
        self.reference_genome = reference_genome

    def create_sequence_dictionary(self, output_dict, picard_path="gatk"):
        """
        Create a sequence dictionary using GATK (Picard runs through GATK)
        """
        if os.path.exists(output_dict):
            print(f"Skipping sequence dictionary creation; file already exists: {output_dict}")
            return
        cmd = [
            picard_path,
            "CreateSequenceDictionary",
            "-R", self.reference_genome,
            "-O", output_dict
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"Sequence dictionary created in file {output_dict}")

    def ensure_fasta_index(self):
        """
        Check for the existence of a fasta index file, and create it if doesn't exist
        """
        fai_path = self.reference_genome + ".fai"
        if not os.path.exists(fai_path):
            print(f"Fasta Index file not found: {fai_path}. Creating index using samtools faidx...")
            cmd = [
                "samtools",
                "faidx",
                self.reference_genome
            ]
            subprocess.run(cmd, check=True)
            print(f"Created Fasta index: {fai_path}")

    def convert_gff_to_interval_list(self, gff_file, dict_file, output_file, gene_positions_file, feature_type='gene', attribute='ID'):
        """
        Convert GFF3 file to a GATK-compatible interval list using gffutils
        """
        header_lines = []
        with open(dict_file, 'r') as f:
            for line in f:
                if line.startswith("@"):
                    header_lines.append(line.strip())
                else:
                    break
        
        db_path = f"{gff_file}.db"
        if not os.path.exists(db_path):
            gffutils.create_db(
                gff_file,
                dbfn=db_path,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True
            )
        db = gffutils.FeatureDB(db_path)

        genes = []
        for feature in db.features_of_type(feature_type):
            gene_id = feature.attributes.get("ID", [None])[0]
            gene_name = feature.attributes.get("gene", [None])[0]
            if gene_name is None:
                gene_name = feature.attributes.get("Name", [gene_id])[0]
            if gene_id:
                genes.append((feature.chrom, feature.start, feature.end, feature.strand, gene_id, gene_name))

        with open(output_file, 'w') as out:
            for header in header_lines:
                out.write(f"{header}\n")
            for chrom, start, end, strand, gene_id, gene_name in genes:
                out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{gene_name}\n")

        with open(gene_positions_file, "w") as gpos:
            for chrom, start, end, strand, gene_id, gene_name in genes:
                gpos.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t{gene_id}\n")
        
        return len(genes)
    
    def validate_bam_file(self, bam_path, gatk_path):
        """
        Perform basic quality control on a BAM file and 
        create BAM index (.bai) using GATK if missing

        Checks:
            - File exists and is readable
            - BAM is indexed (.bai file exists)
            - Contains at least one alignment

        Args:
            bam_path (str): Path to BAM file
            gatk_path (str): Path to GATK executable

        Raises:
            FileNotFoundError: If BAM file doesn't exist
            RuntimeError: If file is unreadable or lacks alignments (empty)
        """
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        # Define expected bam index file path    
        bai_path = f"{bam_path}.bai"
        alt_bai_path = bam_path.replace(".bam", ".bai")

        # If index is missing 
        if not os.path.exists(bai_path) and not os.path.exists(alt_bai_path):
            print(f"BAM index (.bai) file not found for: {bam_path}. Creating index using GATK...")
            cmd = [
                gatk_path,
                "BuildBamIndex",
                "-I", bam_path,
                "--OUTPUT", bai_path
            ]
            subprocess.run(cmd, check=True)
            print(f"Index created for {bam_path}")
        
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            first_read = next(bam.fetch(until_eof=True), None)
            if first_read is None:
                raise ValueError("No alignments found in BAM file.")
            print(f"BAM file '{bam_path}' is valid and contains alignments.")
            bam.close()
        except Exception as e:
            raise RuntimeError(f"Failed to read BAM file: {e}")
        
    def run_depth_of_coverage(self, bam_path, reference_path, interval_list, output_prefix, gatk_path="gatk"):
        """
        Run GATK DepthOfCoverage to generate generate coverage data

        Args:
            bam_path (str): Path to input BAM file
            reference_path (str): Path to reference FASTA
            interval_list (str): Path to interval list file
            output_prefix (str): Output file prefix (without file extension)
            gatk_path (str): Path to GATK executable
        """
        self.ensure_fasta_index()
        cmd = [
            gatk_path,
            "DepthOfCoverage",
            "-R", reference_path,
            "-I", bam_path,
            "-L", interval_list,
            "-O", output_prefix,
            "--omit-depth-output-at-each-base",
            "--omit-locus-table",
            "--summary-coverage-threshold", "10"
        ]
        print("Running GATK DepthOfCoverage...")
        subprocess.run(cmd, check=True)
        print(f"DepthOfCoverage complete. Summary file: {output_prefix}.sample_interval_summary")
