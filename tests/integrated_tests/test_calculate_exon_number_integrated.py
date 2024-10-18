import polars as pl
import pytest
import RNApysoforms as RNApy
import os

def test_integrated_exon_calculation():
    # Get the directory of the current file
    current_dir = os.path.dirname(__file__)

    # Construct the path to the test_data directory
    test_data_dir = os.path.abspath(os.path.join(current_dir, '../test_data/'))

    # Define file paths
    gtf_path = os.path.join(test_data_dir, "Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

    # Read gtf
    annotation = RNApy.read_ensembl_gtf(gtf_path)

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
    try:
        assert sorted_annotation["exon_number"].to_list() == sorted_annotation_2["exon_number"].to_list(), "Exon numbers do not match"
        assert sorted_annotation["start"].to_list() == sorted_annotation_2["start"].to_list(), "Start positions do not match"
        assert sorted_annotation["end"].to_list() == sorted_annotation_2["end"].to_list(), "End positions do not match"

        print("Test succeeded: Exon numbers and coordinates match.")

    except AssertionError as e:
        print(f"Test failed: {e}")
        raise  # Re-raise exception to let pytest know the test failed
