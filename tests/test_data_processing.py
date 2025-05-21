import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from gene_copy_number.data_processing import DataPreprocessor

def test_data_processing():
    test_dir = os.path.join(os.path.dirname(__file__), "../test_data")
    reference = "fake_reference.fasta"
    gff = os.path.join(test_dir, "test_annotations.gff3")
    dict_file = os.path.join(test_dir, "reference.dict")
    output_interval_list = os.path.join(test_dir, "test_intervals.list")
    bam_file = os.path.join(test_dir, "test.bam")

    try:
        dp = DataPreprocessor(reference)
        dp.validate_bam_file(bam_file)  # will auto-create index if missing
        count = dp.convert_gff_to_interval_list(
            gff_file=gff,
            dict_file=dict_file,
            output_file=output_interval_list
        )
        assert count == 2, "Expected 2 genes in interval list"
        assert os.path.exists(output_interval_list), "Interval list not generated"
        print("✅ Data processing test passed.")
    except Exception as e:
        print(f"❌ Data processing test failed: {e}")

if __name__ == "__main__":
    test_data_processing()