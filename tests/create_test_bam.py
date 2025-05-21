import pysam
import os

bam_path = "test_data/test.bam"
os.makedirs("test_data", exist_ok=True)

header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}
with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
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

pysam.index(bam_path)
print(f"✅ Created BAM file and index at {bam_path}")