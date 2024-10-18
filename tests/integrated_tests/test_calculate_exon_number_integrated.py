import polars as pl
import pytest
import RNApysoforms as pt
import os


def integrated_exon_calculation_test():

    # Get the directory of the current file
    current_dir = os.path.dirname(__file__)

    # Construct the path to the test_data directory
    test_data_dir = os.path.abspath(os.path.join(current_dir, '../test_data/'))

    # Define file paths
    gtf_path = os.path.join(test_data_dir, "Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

    # Read gtf
    annotation = pt.read_ensembl_gtf("./test_data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

    # Add introns
    annotation = pt.to_intron(annotation)

    # Drop exon number
    annotation_2 = pt.annotation.drop("exon_number")

    # Recalculate exon number
    annotation_2 = pt.calculate_exon_number(annotation_2)

    ## Check if recalculated exon numbers and coordinates reflect the original.
    assert annotation.sort(["transcript_id", "start"])["exon_number"].to_list() == \
    annotation_2.sort(["transcript_id", "start"])["exon_number"].to_list()

    assert annotation.sort(["transcript_id", "start"])["start"].to_list() == \
        annotation_2.sort(["transcript_id", "start"])["start"].to_list()

    assert annotation.sort(["transcript_id", "start"])["end"].to_list() == \
        annotation_2.sort(["transcript_id", "start"])["end"].to_list()
