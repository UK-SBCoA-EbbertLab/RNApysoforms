# tests/test_read_ensembl_gtf.py

import pytest
import polars as pl
import tempfile
import os
from RNApysoforms import read_ensembl_gtf

def test_read_ensembl_gtf_valid_file():
    """
    Test reading a valid GTF file with expected content.
    """
    # Create a sample GTF content
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; transcript_name "DDX11L1-002"; transcript_biotype "transcribed_unprocessed_pseudogene"; exon_number "1";
chr1\tHAVANA\tCDS\t12010\t12057\t.\t+\t0\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; transcript_name "DDX11L1-002"; transcript_biotype "transcribed_unprocessed_pseudogene"; exon_number "1";
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; transcript_name "DDX11L1-002"; transcript_biotype "transcribed_unprocessed_pseudogene";
"""
    # Write the content to a temporary GTF file
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        # Call the function
        df = read_ensembl_gtf(tmp_gtf_path)

        # Verify that the returned object is a Polars DataFrame
        assert isinstance(df, pl.DataFrame)

        # Verify that the DataFrame has the expected columns
        expected_columns = [
            "gene_id",
            "gene_name",
            "transcript_id",
            "transcript_name",
            "transcript_biotype",
            "seqnames",
            "strand",
            "type",
            "start",
            "end",
            "exon_number"
        ]
        assert df.columns == expected_columns

        # Verify that the DataFrame has the expected number of rows (only 'exon' and 'CDS' features)
        assert df.shape[0] == 2  # Only the first two entries are 'exon' and 'CDS'

        # Verify the content of the DataFrame
        expected_data = {
            "gene_id": ["ENSG00000223972", "ENSG00000223972"],
            "gene_name": ["DDX11L1", "DDX11L1"],
            "transcript_id": ["ENST00000456328", "ENST00000456328"],
            "transcript_name": ["DDX11L1-002", "DDX11L1-002"],
            "transcript_biotype": ["transcribed_unprocessed_pseudogene", "transcribed_unprocessed_pseudogene"],
            "seqnames": ["chr1", "chr1"],
            "strand": ["+", "+"],
            "type": ["exon", "CDS"],
            "start": [11869, 12010],
            "end": [12227, 12057],
            "exon_number": [1, 1]
        }
        expected_df = pl.DataFrame(expected_data)
        assert df.equals(expected_df, null_equal=True)

    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_nonexistent_file():
    """
    Test that ValueError is raised when the file path does not exist.
    """
    with pytest.raises(ValueError) as excinfo:
        read_ensembl_gtf("nonexistent_file.gtf")
    assert "File 'nonexistent_file.gtf' does not exist." in str(excinfo.value)

def test_read_ensembl_gtf_path_is_directory():
    """
    Test that ValueError is raised when the path is not a file.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        with pytest.raises(ValueError) as excinfo:
            read_ensembl_gtf(tmp_dir)
        assert f"'{tmp_dir}' is not a file." in str(excinfo.value)

def test_read_ensembl_gtf_invalid_extension():
    """
    Test that ValueError is raised when the file does not have a '.gtf' extension.
    """
    # Create a temporary file with an invalid extension
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as tmp_file:
        tmp_file.write("Sample content")
        tmp_file_path = tmp_file.name

    try:
        with pytest.raises(ValueError) as excinfo:
            read_ensembl_gtf(tmp_file_path)
        assert "File must have a '.gtf' extension." in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(tmp_file_path)

def test_read_ensembl_gtf_missing_columns():
    """
    Test that ValueError is raised when required columns are missing or file cannot be read properly.
    """
    # Create a GTF file with missing columns
    gtf_content = "Invalid content without proper columns"
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        with pytest.raises(Exception) as excinfo:
            read_ensembl_gtf(tmp_gtf_path)

        print(excinfo.value)
        # The error message can vary; checking for general failure
        assert "The length of the new names list should be equal to or less than the original column length" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_no_exon_or_cds_features():
    """
    Test that the function returns an empty DataFrame when there are no 'exon' or 'CDS' features.
    """
    # Create a GTF file without 'exon' or 'CDS' features
    gtf_content = """\
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that the DataFrame is empty
        assert df.is_empty()
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_missing_gene_name_and_transcript_name():
    """
    Test that 'gene_name' and 'transcript_name' are filled with 'gene_id' and 'transcript_id' when missing.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that 'gene_name' and 'transcript_name' are filled
        assert df['gene_name'][0] == "ENSG00000223972"
        assert df['transcript_name'][0] == "ENST00000456328"
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_missing_exon_number():
    """
    Test that 'exon_number' is handled correctly when missing.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; transcript_name "DDX11L1-002";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that 'exon_number' is null or NaN
        assert df['exon_number'][0] is None or pl.is_nan(df['exon_number'][0])
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_null_values_after_parsing():
    """
    Test that ValueError is raised when there are null values after parsing attributes.
    """
    # Create a GTF file with missing 'gene_id' in attributes
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\ttranscript_id "ENST00000456328"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        with pytest.raises(ValueError) as excinfo:
            read_ensembl_gtf(tmp_gtf_path)
        assert "This GTF file is not consistent with the 2024 ENSEMBL GTF format." in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_additional_attributes():
    """
    Test with a GTF file containing additional or unexpected attributes.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "1"; gene_name "DDX11L1"; extra_attribute "value"; transcript_id "ENST00000456328"; transcript_name "DDX11L1-002"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that standard attributes are parsed correctly
        assert df['gene_id'][0] == "ENSG00000223972"
        assert df['gene_name'][0] == "DDX11L1"
        # Verify that extra attributes are ignored
        assert 'extra_attribute' not in df.columns
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_attributes_different_order():
    """
    Test with a GTF file where attributes are in a different order.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\ttranscript_name "DDX11L1-002"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; gene_id "ENSG00000223972"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that attributes are parsed correctly despite the order
        assert df['gene_id'][0] == "ENSG00000223972"
        assert df['gene_name'][0] == "DDX11L1"
        assert df['transcript_id'][0] == "ENST00000456328"
        assert df['transcript_name'][0] == "DDX11L1-002"
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_correct_columns_and_types():
    """
    Test that the DataFrame has the correct columns and data types.
    """
    gtf_content = """\
chr1\tHAVANA\tCDS\t12010\t12057\t.\t+\t0\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify columns
        expected_columns = [
            "gene_id",
            "gene_name",
            "transcript_id",
            "transcript_name",
            "transcript_biotype",
            "seqnames",
            "strand",
            "type",
            "start",
            "end",
            "exon_number"
        ]
        assert df.columns == expected_columns

        # Verify data types
        expected_dtypes = {
            "gene_id": pl.Utf8,
            "gene_name": pl.Utf8,
            "transcript_id": pl.Utf8,
            "transcript_name": pl.Utf8,
            "transcript_biotype": pl.Utf8,
            "seqnames": pl.Utf8,
            "strand": pl.Utf8,
            "type": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "exon_number": pl.Int64
        }
        for col in df.columns:
            assert df.schema[col] == expected_dtypes[col]
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_filters_out_other_features():
    """
    Test that the function filters out features other than 'exon' and 'CDS'.
    """
    gtf_content = """\
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\ttranscript_id "ENST00000456328"; gene_id "ENSG00000223972";
chr1\tHAVANA\tstart_codon\t12010\t12012\t.\t+\t0\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chr1\tHAVANA\tCDS\t12010\t12057\t.\t+\t0\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Only 'CDS' feature should be included
        assert df.shape[0] == 1
        assert df['type'][0] == 'CDS'
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_handles_phase_and_score():
    """
    Test that the function handles different 'phase' and 'score' values.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t1000\t+\t2\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that 'score' and 'phase' columns are handled (though not included in the final DataFrame)
        assert df['start'][0] == 11869
        assert df['end'][0] == 12227
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)

def test_read_ensembl_gtf_handles_strand_values():
    """
    Test that the function handles different 'strand' values, including invalid ones.
    """
    gtf_content = """\
chr1\tHAVANA\texon\t11869\t12227\t.\t*\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
chr1\tHAVANA\texon\t12500\t12700\t.\t-\t.\tgene_id "ENSG00000223973"; transcript_id "ENST00000456329"; exon_number "1";
"""
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gtf', delete=False) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        tmp_gtf_path = tmp_gtf.name

    try:
        df = read_ensembl_gtf(tmp_gtf_path)
        # Verify that both entries are included
        assert df.shape[0] == 2
        # Check 'strand' values
        assert df['strand'][0] == '*'
        assert df['strand'][1] == '-'
    finally:
        # Clean up the temporary file
        os.remove(tmp_gtf_path)
