import polars as pl
import pytest
import RNApysoforms as RNApy
import os

@pytest.fixture
def get_test_data():
    current_dir = os.path.dirname(__file__)
    test_data_dir = os.path.abspath(os.path.join(current_dir, '../test_data/'))
    gtf_path = os.path.join(test_data_dir, "Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")
    return gtf_path

def test_integrated_exon_calculation(get_test_data):
    # Read gtf
    annotation = RNApy.read_ensembl_gtf(get_test_data)

    # Add introns
    annotation = RNApy.to_intron(annotation)

    # Drop exon number
    annotation_2 = annotation.drop("exon_number")

    # Recalculate exon number
    annotation_2 = RNApy.calculate_exon_number(annotation_2)

    # Sort annotations
    sorted_annotation = annotation.sort(["transcript_id", "start"])
    sorted_annotation_2 = annotation_2.sort(["transcript_id", "start"])

    # Check if recalculated exon numbers and coordinates reflect the original.
    assert sorted_annotation["exon_number"].to_list() == sorted_annotation_2["exon_number"].to_list()
    assert sorted_annotation["start"].to_list() == sorted_annotation_2["start"].to_list()
    assert sorted_annotation["end"].to_list() == sorted_annotation_2["end"].to_list()

    # Additional assertions
    # Check other columns like 'strand', 'gene_id', etc.
    assert sorted_annotation["strand"].to_list() == sorted_annotation_2["strand"].to_list()
    assert sorted_annotation["gene_id"].to_list() == sorted_annotation_2["gene_id"].to_list()

@pytest.mark.parametrize("test_input", [
    "single_exon_transcript.gtf",
    "overlapping_exons.gtf",
    "empty_annotation.gtf",
    # Add paths to other test GTF files
])
def test_edge_cases(test_input):
    # Assuming test_input is the path to a test GTF file
    annotation = RNApy.read_ensembl_gtf(test_input)

    # Proceed with the same steps as before
    # ...

    # Include try-except blocks to test error handling
    try:
        # Attempt to process annotation
        pass
    except Exception as e:
        assert isinstance(e, ExpectedExceptionType)