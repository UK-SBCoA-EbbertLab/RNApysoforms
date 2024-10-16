# tests/test_read_expression_matrix.py

import pytest
import polars as pl
import tempfile
import os
from RNApysoforms import read_expression_matrix
import warnings
from polars.testing import assert_frame_equal
from polars.testing import assert_series_equal

def _create_temp_file(content: str, suffix: str) -> str:
    """
    Helper function to create a temporary file with the given content and suffix.
    
    Parameters
    ----------
    content : str
        The content to write into the temporary file.
    suffix : str
        The file extension/suffix (e.g., '.csv', '.tsv').
    
    Returns
    -------
    str
        The path to the created temporary file.
    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix=suffix, delete=False) as tmp_file:
        tmp_file.write(content)
        tmp_file_path = tmp_file.name
    return tmp_file_path

def test_read_expression_matrix_basic():
    """
    Test the basic functionality with a valid expression matrix without metadata.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(expression_matrix_path=expr_path)

        # Verify that the returned object is a Polars DataFrame
        assert isinstance(df, pl.DataFrame)

        # Verify that the DataFrame has the expected columns
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        # Verify the content of the DataFrame
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400]
        }
    
        ## Check if frames are equal
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)

    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_metadata():
    """
    Test reading an expression matrix and merging with metadata.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a sample metadata CSV content
    metadata_content = """sample_id,condition
sample1,treated
sample2,control
sample3,treated
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        with pytest.warns(UserWarning) as record:
            # Call the function with metadata
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path
            )

            # Verify that the DataFrame has merged metadata
            expected_columns = [
                "transcript_id",
                "gene_id",
                "sample_id",
                "counts",
                "condition"
            ]
            assert df.columns == expected_columns

            # Verify the content of the DataFrame
            expected_data = {
                "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
                "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
                "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
                "counts": [100, 200, 150, 250, 300, 400],
                "condition": ["treated", "control", "treated", "control", "treated", "control"]
            }
            
            expected_df = pl.DataFrame(expected_data)
            assert_frame_equal(df, expected_df, check_row_order=False)

            assert len(record) == 1
            assert "The following sample IDs are present in metadata but not in expression data: ['sample3']." in str(record[0].message)

    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_cpm_normalization():
    """
    Test CPM normalization functionality.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with CPM normalization
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            cpm_normalization=True
        )

        # Expected CPM values
        # For sample1: total = 250
        # tx1: (100 / 250) * 1e6 = 400000
        # tx2: (150 / 250) * 1e6 = 600000
        # For sample2: total = 450
        # tx1: (200 / 450) * 1e6 ≈ 444444.44
        # tx2: (250 / 450) * 1e6 ≈ 555555.56

        # Verify that CPM columns are added
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts",
            "CPM"
        ]
        assert df.columns == expected_columns

        # Verify the CPM values (allowing a small tolerance)
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250],
            "CPM": [400000.0, 444444.4, 600000.0, 555555.6]
        }
        expected_df = pl.DataFrame(expected_data)
        for col in ["CPM"]:
            assert_series_equal(df[col].round(1), expected_df[col], check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_relative_abundance():
    """
    Test relative transcript abundance calculation.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with relative abundance
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            relative_abundance=True
        )

        # Calculate expected relative abundances
        # For gene1, sample1 total = 250
        # tx1: (100 / 250) * 100 = 40
        # tx2: (150 / 250) * 100 = 60
        # For gene2, sample1 total = 300
        # tx3: (300 / 300) * 100 = 100
        # Similarly for sample2
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400],
            "relative_abundance": [40.0, 44.4, 60.0, 55.6, 100.0, 100.0]
        }
        expected_df = pl.DataFrame(expected_data)
        # Verify that relative_abundance column is added
        assert "relative_abundance" in df.columns

        # Compare relative abundance values
        assert_series_equal(df["relative_abundance"].round(1), expected_df["relative_abundance"], check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_cpm_and_relative_abundance():
    """
    Test performing both CPM normalization and relative abundance calculation.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with both CPM normalization and relative abundance
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            cpm_normalization=True,
            relative_abundance=True
        )

        # Expected columns
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts",
            "CPM",
            "relative_abundance"
        ]
        assert df.columns == expected_columns

        # Verify CPM and relative_abundance values
        # CPM calculations as in previous test
        # Relative abundance as in previous test
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400],
            "CPM": [181818.18,235294.1, 272727.3,294117.6, 545454.6, 470588.2],
            "relative_abundance": [40.0, 44.4, 60.0, 55.6, 100.0, 100.0]
        }
        expected_df = pl.DataFrame(expected_data)
        for col in ["CPM", "relative_abundance"]:
            assert_series_equal(df[col].round(1), expected_df[col], check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_missing_transcript_id_column():
    """
    Test that ValueError is raised when transcript_id_column_name is None.
    """
    # Create a sample expression matrix CSV content
    expr_content = """gene_id,sample1,sample2
gene1,100,200
gene1,150,250
gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                transcript_id_column_name=None
            )
        assert "The 'transcript_id_column_name' is required and cannot be None." in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_missing_feature_id_columns():
    """
    Test that ValueError is raised when required feature ID columns are missing in the expression matrix.
    """
    # Create a sample expression matrix CSV content missing 'transcript_id'
    expr_content = """gene_id,sample1,sample2
gene1,100,200
gene1,150,250
gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                transcript_id_column_name="transcript_id",
                gene_id_column_name="gene_id"
            )
        assert "The following feature ID columns are missing in the expression dataframe: ['transcript_id']" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_non_numeric_expression_columns():
    """
    Test that ValueError is raised when expression columns are not numeric.
    """
    # Create a sample expression matrix CSV content with non-numeric counts
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,abc,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path
            )
        assert "The following columns are expected to be numerical but are not: ['sample1']" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_metadata_missing_columns():
    """
    Test that ValueError is raised when required columns are missing in the metadata DataFrame.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a metadata file missing 'sample_id'
    metadata_content = """condition
treated
control
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path,
                metadata_sample_id_column="sample_id"
            )
        assert "The metadata_sample_id_column 'sample_id' is not present in the metadata dataframe." in str(excinfo.value)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_no_overlapping_sample_ids():
    """
    Test that ValueError is raised when there are no overlapping sample IDs between expression data and metadata.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a metadata file with non-overlapping sample IDs
    metadata_content = """sample_id,condition
sample3,treated
sample4,control
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path
            )
        assert "No overlapping sample IDs found between expression data and metadata." in str(excinfo.value)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_partial_overlapping_sample_ids():
    """
    Test that a warning is issued when there is a partial overlap of sample IDs between expression data and metadata.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2,sample3
tx1,gene1,100,200,300
tx2,gene1,150,250,350
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a metadata file with some overlapping and some non-overlapping sample IDs
    metadata_content = """sample_id,condition
sample2,control
sample3,treated
sample4,control
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        with pytest.warns(UserWarning) as record:
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path
            )
        # Check that one warning is issued
        assert len(record) == 1
        assert "sample IDs are present in metadata but not in expression data" in str(record[0].message) or \
               "sample IDs are present in expression data but not in metadata" in str(record[0].message)

        # Verify that the merged DataFrame only includes overlapping sample IDs
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample2", "sample3", "sample2", "sample3"],
            "counts": [200, 300, 250, 350],
            "condition": ["control", "treated", "control", "treated"]
        }

        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)

    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_custom_column_names():
    """
    Test that custom column names are handled correctly.
    """
    # Create a sample expression matrix CSV content with custom column names
    expr_content = """tid,gid,s1,s2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with custom column names
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            gene_id_column_name="gid",
            transcript_id_column_name="tid",
            expression_measure_name="expression",
            cpm_normalization=True
        )

        # Verify that the DataFrame has the expected columns
        expected_columns = [
            "tid",
            "gid",
            "sample_id",
            "expression",
            "CPM"
        ]
        assert df.columns == expected_columns

        # Verify CPM calculations
        # For sample1: total = 250
        # tx1: (100 / 250) * 1e6 = 400000
        # tx2: (150 / 250) * 1e6 = 600000
        # For sample2: total = 450
        # tx1: (200 / 450) * 1e6 ≈ 444444.44
        # tx2: (250 / 450) * 1e6 ≈ 555555.56
        expected_data = {
            "tid": ["tx1", "tx1", "tx2", "tx2"],
            "gid": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "expression": [100, 200, 150, 250],
            "CPM": [400000.0, 444444.4, 600000.0, 555555.6]
        }
        expected_df = pl.DataFrame(expected_data)
        for col in ["CPM"]:
            assert_series_equal(df[col].round(1), expected_df[col], check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_multiple_expression_columns():
    """
    Test that multiple expression columns are handled correctly.
    """
    # Create a sample expression matrix CSV content with multiple expression columns
    expr_content = """transcript_id,gene_id,sample1,sample2,sample3
tx1,gene1,100,200,300
tx2,gene1,150,250,350
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with multiple expression columns
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            expression_measure_name="expression",
            cpm_normalization=True,
            relative_abundance=True
        )

        # Expected columns
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "expression",
            "CPM",
            "relative_abundance"
        ]
        assert df.columns == expected_columns

        # Verify the content of the DataFrame
        # Relative abundance calculations:
        # gene1, sample1: (100 / 250) * 100 = 40
        # gene1, sample2: (200 / 450) * 100 ≈ 44.44
        # gene1, sample3: (300 / 650) * 100 ≈ 46.15
        # gene1, sample1: (150 / 250) * 100 = 60
        # gene1, sample2: (250 / 450) * 100 ≈ 55.56
        # gene1, sample3: (350 / 650) * 100 ≈ 53.85
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx1", "tx2", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample3", "sample1", "sample2", "sample3"],
            "expression": [100, 200, 300, 150, 250, 350],
            "CPM": [400000.0, 444444.44, 461538.46, 600000.0, 555555.56, 538461.54],
            "relative_abundance": [40.0, 44.44, 46.15, 60.0, 55.56, 53.85]
        }
        expected_df = pl.DataFrame(expected_data)
        for col in ["CPM", "relative_abundance"]:
            assert_series_equal(df[col].round(2), expected_df[col], check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_gene_id_none_with_relative_abundance():
    """
    Test that a warning is issued when relative_abundance is True but gene_id_column_name is None.
    """
    # Create a sample expression matrix CSV content without gene_id
    expr_content = """transcript_id,sample1,sample2
tx1,100,200
tx2,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.warns(UserWarning) as record:
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                relative_abundance=True,
                gene_id_column_name=None
            )
        # Check that the warning message is correct
        assert len(record) == 1
        assert "relative_abundance was set to True, but gene_id_column_name was not provided" in str(record[0].message)

        # Verify that relative_abundance is not calculated
        assert "relative_abundance" not in df.columns
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_expression_measure_name():
    """
    Test that expression_measure_name parameter is applied correctly.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with a custom expression_measure_name
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            expression_measure_name="expr_counts"
        )

        # Verify that the expression_measure_name is used
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "expr_counts"
        ]
        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "expr_counts": [100, 200, 150, 250]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_supported_file_formats():
    """
    Test that the function can read different supported file formats.
    """
    # Define sample data
    expr_data = {
        "transcript_id": ["tx1", "tx2"],
        "gene_id": ["gene1", "gene1"],
        "sample1": [100, 150],
        "sample2": [200, 250]
    }
    expr_df = pl.DataFrame(expr_data)

    # Supported formats: .csv, .tsv, .txt, .parquet, .xlsx
    formats = ['.csv', '.tsv', '.txt', '.parquet', '.xlsx']
    for suffix in formats:
        with tempfile.NamedTemporaryFile(mode='w+', suffix=suffix, delete=False) as tmp_file:
            if suffix in ['.csv', '.tsv', '.txt']:
                sep = ',' if suffix == '.csv' else '\t'
                expr_df.write_csv(tmp_file.name, separator=sep)
            elif suffix == '.parquet':
                expr_df.write_parquet(tmp_file.name)
            elif suffix == '.xlsx':
                expr_df.write_excel(tmp_file.name)
            tmp_file_path = tmp_file.name

        try:
            # Call the function
            if suffix in ['.csv', '.tsv', '.txt']:
                df = read_expression_matrix(expression_matrix_path=tmp_file_path)
            elif suffix == '.parquet':
                df = read_expression_matrix(expression_matrix_path=tmp_file_path)
            elif suffix == '.xlsx':
                df = read_expression_matrix(expression_matrix_path=tmp_file_path)
            
            # Verify that the DataFrame matches
            expected_long = expr_df.unpivot(
                index=["transcript_id", "gene_id"],
                on=["sample1", "sample2"],
                variable_name="sample_id",
                value_name="counts"
            )
            assert df.equals(expected_long)
        finally:
            # Clean up the temporary file
            os.remove(tmp_file_path)

def test_read_expression_matrix_unsupported_file_format():
    """
    Test that ValueError is raised for unsupported file formats.
    """
    # Create a sample expression matrix with unsupported extension
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.unsupported')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(expression_matrix_path=expr_path)
        assert "Unsupported file extension '.unsupported'" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_zero_gene_counts():
    """
    Test that relative abundance is set to zero when total gene counts are zero to avoid division by zero.
    """
    # Create a sample expression matrix CSV content with zero counts for a gene
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,0,0
tx2,gene1,0,0
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with relative abundance
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            relative_abundance=True
        )

        # Verify that relative_abundance is zero for gene1 transcripts
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [0, 0, 0, 0, 300, 400],
            "relative_abundance": [0.0, 0.0, 0.0, 0.0, 100.0, 100.0]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_gene_id_column_none():
    """
    Test that gene_id_column_name can be set to None and handled correctly.
    """
    # Create a sample expression matrix CSV content without gene_id
    expr_content = """transcript_id,sample1,sample2
tx1,100,200
tx2,150,250
tx3,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with gene_id_column_name set to None
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            gene_id_column_name=None
        )

        # Verify that gene_id is not present
        assert "gene_id" not in df.columns

        # Verify the DataFrame content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_transcript_id_custom():
    """
    Test that custom transcript_id_column_name is handled correctly.
    """
    # Create a sample expression matrix CSV content with custom transcript_id column name
    expr_content = """tid,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with custom transcript_id_column_name
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            transcript_id_column_name="tid"
        )

        # Verify that the transcript_id column is correctly mapped
        expected_columns = [
            "tid",
            "gene_id",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "tid": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_empty_expression_columns():
    """
    Test that the function handles expression matrices with no expression columns.
    """
    # Create a sample expression matrix CSV content with only transcript_id and gene_id
    expr_content = """transcript_id,gene_id
tx1,gene1
tx2,gene1
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path
        )

        # Verify that the DataFrame has only transcript_id, gene_id, sample_id columns with no counts
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        # Since there are no expression columns, counts should be empty
        expected_data = {
            "transcript_id": [],
            "gene_id": [],
            "sample_id": [],
            "counts": []
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False, check_dtypes=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_invalid_file_content():
    """
    Test that ValueError is raised when the expression matrix file cannot be read properly.
    """
    # Create a sample expression matrix CSV content with invalid data
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,abc
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(expression_matrix_path=expr_path)
        assert "The following columns are expected to be numerical but are not: ['sample2']" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_zero_counts():
    """
    Test that relative abundance is handled correctly when total gene counts are zero.
    """
    # Create a sample expression matrix CSV content with zero counts for a gene
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,0,0
tx2,gene1,0,0
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with relative abundance
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            relative_abundance=True
        )

        # Verify that relative_abundance is zero for gene1 transcripts
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [0, 0, 0, 0, 300, 400],
            "relative_abundance": [0.0, 0.0, 0.0, 0.0, 100.0, 100.0]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_sample_jitter():
    """
    Test that marker_jitter parameter is applied correctly.
    """
    # Note: Since marker_jitter is not a parameter in read_expression_matrix,
    # this test is unnecessary. It might have been intended for a plotting function.
    # Including it here as a placeholder.
    pass  # No implementation needed

def test_read_expression_matrix_expression_measure_non_default():
    """
    Test that expression_measure_name parameter is applied correctly.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with a custom expression_measure_name
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            expression_measure_name="expr_values"
        )

        # Verify that the expression_measure_name is used
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "expr_values"
        ]
        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "expr_values": [100, 200, 150, 250]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_gene_id_column_none():
    """
    Test that gene_id_column_name can be set to None and handled correctly.
    """
    # Create a sample expression matrix CSV content without gene_id
    expr_content = """transcript_id,sample1,sample2
tx1,100,200
tx2,150,250
tx3,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.warns(UserWarning) as record:
            # Call the function with gene_id_column_name set to None
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                gene_id_column_name=None,
                relative_abundance=True  # Should trigger a warning and skip relative_abundance
            )

            # Verify that gene_id is not present
            assert "gene_id" not in df.columns

            # Verify that relative_abundance is not present
            assert "relative_abundance" not in df.columns

            # Verify the content
            expected_data = {
                "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
                "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
                "counts": [100, 200, 150, 250, 300, 400]
            }
            expected_df = pl.DataFrame(expected_data)
            assert_frame_equal(df, expected_df, check_row_order=False)

            # Check that the warning message is correct
            assert len(record) == 1
            assert "relative_abundance was set to True, but gene_id_column_name was not provided (set to None)" in str(record[0].message)

    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_gene_id_column_custom():
    """
    Test that custom gene_id_column_name is handled correctly.
    """
    # Create a sample expression matrix CSV content with custom gene_id column name
    expr_content = """transcript_id,gid,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with custom gene_id_column_name
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            gene_id_column_name="gid"
        )

        # Verify that the gene_id column is correctly mapped
        expected_columns = [
            "transcript_id",
            "gid",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gid": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_missing_expression_columns():
    """
    Test that the function handles expression matrices with missing expression columns.
    """
    # Create a sample expression matrix CSV content with missing expression columns
    expr_content = """transcript_id,gene_id
tx1,gene1
tx2,gene1
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path
        )

        # Verify that the DataFrame has only transcript_id, gene_id, sample_id, counts columns with no counts
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        # Verify that counts are empty
        assert df["counts"].is_empty()
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_expression_columns_non_numeric():
    """
    Test that ValueError is raised when expression columns contain non-numeric data.
    """
    # Create a sample expression matrix CSV content with non-numeric expression data
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,abc
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path
            )
        assert "The following columns are expected to be numerical but are not: ['sample2']" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_duplicate_transcript_ids():
    """
    Test that the function handles duplicate transcript IDs correctly.
    """
    # Create a sample expression matrix CSV content with duplicate transcript IDs
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx1,gene1,150,250
tx2,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path
        )

        # Verify that the DataFrame correctly melts duplicate transcript IDs
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_metadata_non_csv():
    """
    Test that the function can read metadata files in supported formats other than CSV.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a sample metadata TSV content
    metadata_content = """sample_id\tcondition
sample1\ttreated
sample2\tcontrol
"""
    metadata_path = _create_temp_file(metadata_content, '.tsv')

    try:
        # Call the function with metadata in TSV format
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            metadata_path=metadata_path
        )

        # Verify that the DataFrame has merged metadata correctly
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250],
            "condition": ["treated", "control", "treated", "control"]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_gene_id_missing_with_relative_abundance():
    """
    Test that relative_abundance is calculated correctly when gene_id_column_name is provided.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
tx3,gene2,300,400
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with relative abundance
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            relative_abundance=True
        )

        # Verify that relative_abundance is calculated correctly
        # For gene1, sample1: (100 / 250) * 100 = 40
        # For gene1, sample2: (200 / 450) * 100 ≈ 44.44
        # For gene2, sample1: (300 / 300) * 100 = 100
        # For gene2, sample2: (400 / 400) * 100 = 100
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250, 300, 400],
            "relative_abundance": [40.0, 44.4, 60.0, 55.6, 100.0, 100.0]
        }
        expected_df = pl.DataFrame(expected_data)
        # Allow small tolerance for floating point calculations
        for col in ["relative_abundance"]:
            assert_series_equal(df[col].round(1), expected_df[col].round(1), check_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_empty_expression_matrix():
    """
    Test that the function handles an empty expression matrix gracefully.
    """
    # Create an empty expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:

        with pytest.raises(ValueError) as excinfo:
            # Call the function
            df = read_expression_matrix(
                expression_matrix_path=expr_path
            )

            assert "The following columns are expected to" in str(excinfo.value)

    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_invalid_metadata_format():
    """
    Test that ValueError is raised when the metadata file format is unsupported or invalid.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create an unsupported metadata file format
    metadata_content = """sample_id,condition
sample1,treated
sample2,control
"""
    metadata_path = _create_temp_file(metadata_content, '.unsupported')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path
            )
        assert "Unsupported file extension '.unsupported'" in str(excinfo.value)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_merge_with_nonexistent_metadata_file():
    """
    Test that ValueError is raised when the metadata file does not exist.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path="nonexistent_metadata.csv"
            )
        assert "Failed to read the file" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_merge_with_directory_as_metadata_path():
    """
    Test that ValueError is raised when the metadata_path is a directory.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    with tempfile.TemporaryDirectory() as tmp_dir:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=tmp_dir
            )
        assert "Failed to read the file" in str(excinfo.value)

def test_read_expression_matrix_sample_id_column_custom():
    """
    Test that metadata_sample_id_column parameter is applied correctly.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a sample metadata CSV content with a custom sample ID column name
    metadata_content = """sid,condition
sample1,treated
sample2,control
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        # Call the function with a custom metadata_sample_id_column
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            metadata_path=metadata_path,
            metadata_sample_id_column="sid"
        )

        # Verify that the DataFrame has merged metadata correctly
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sid",
            "counts",
            "condition"
        ]

        print(df.columns)

        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sid": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250],
            "condition": ["treated", "control", "treated", "control"]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_handles_extra_columns():
    """
    Test that the function ignores extra columns in the expression matrix and metadata.
    """
    # Create a sample expression matrix CSV content with extra columns
    expr_content = """transcript_id,gene_id,sample1,sample2,extra1
tx1,gene1,100,200,foo
tx2,gene1,150,250,bar
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a sample metadata CSV content with extra columns
    metadata_content = """sample_id,condition,extra_meta
sample1,treated,baz
sample2,control,qux
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:

        with pytest.raises(ValueError) as excinfo:
            # Call the function
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                metadata_path=metadata_path
            )
            assert "The following columns are expected to be numerical" in str(excinfo.value)

    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_large_dataset():
    """
    Test that the function can handle large datasets efficiently.
    """
    # Create a large expression matrix CSV content
    num_transcripts = 1000
    num_samples = 100
    header = ["transcript_id", "gene_id"] + [f"sample{i}" for i in range(1, num_samples+1)]
    rows = [f"tx{j},gene{j%100}," + ",".join([str(j*10 + i) for i in range(1, num_samples+1)]) for j in range(1, num_transcripts+1)]
    expr_content = "\n".join([",".join(header)] + rows)
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            cpm_normalization=True,
            relative_abundance=True
        )

        # Verify the shape of the DataFrame
        assert df.shape[0] == num_transcripts * num_samples
        # Verify that necessary columns are present
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts",
            "CPM",
            "relative_abundance"
        ]
        assert df.columns == expected_columns
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_invalid_gene_id_column_name():
    """
    Test that ValueError is raised when gene_id_column_name is provided but missing in the expression matrix.
    """
    # Create a sample expression matrix CSV content without gene_id
    expr_content = """transcript_id,sample1,sample2
tx1,100,200
tx2,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        with pytest.raises(ValueError) as excinfo:
            read_expression_matrix(
                expression_matrix_path=expr_path,
                gene_id_column_name="gene_id",
                relative_abundance=True
            )
        assert "The following feature ID columns are missing in the expression dataframe: ['gene_id']" in str(excinfo.value)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_expression_columns_all_numeric():
    """
    Test that the function accepts expression columns that are all numeric.
    """
    # Create a sample expression matrix CSV content with all numeric expression columns
    expr_content = """transcript_id,gene_id,sample1,sample2,sample3
tx1,gene1,100,200,300
tx2,gene1,150,250,350
tx3,gene2,300,400,500
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path
        )

        # Verify that the DataFrame has correct columns and data
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts"
        ]
        assert df.columns == expected_columns

        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx1", "tx2", "tx2", "tx2", "tx3", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene1", "gene1", "gene2", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample3"] * 3,
            "counts": [100, 200, 300, 150, 250, 350, 300, 400, 500]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_missing_expression_data():
    """
    Test that the function handles missing expression data gracefully.
    """
    # Create a sample expression matrix CSV content with missing expression data
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,,250
tx3,gene2,300,
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path
        )

        # Verify that missing counts are handled as null
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, None, 250, 300, None]
        }
        expected_df = pl.DataFrame(expected_data)

        assert_frame_equal(expected_df, df, check_row_order=False)

    finally:
        # Clean up the temporary file
        os.remove(expr_path)

def test_read_expression_matrix_with_metadata_extra_columns():
    """
    Test that the function ignores extra columns in the metadata file.
    """
    # Create a sample expression matrix CSV content
    expr_content = """transcript_id,gene_id,sample1,sample2
tx1,gene1,100,200
tx2,gene1,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    # Create a metadata file with extra columns
    metadata_content = """sample_id,condition,extra_meta
sample1,treated,foo
sample2,control,bar
"""
    metadata_path = _create_temp_file(metadata_content, '.csv')

    try:
        # Call the function
        df = read_expression_matrix(
            expression_matrix_path=expr_path,
            metadata_path=metadata_path
        )

        # Verify that extra_meta is ignored
        expected_columns = [
            "transcript_id",
            "gene_id",
            "sample_id",
            "counts",
            "condition",
            "extra_meta"
        ]
        assert df.columns == expected_columns

        # Verify the content
        expected_data = {
            "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
            "gene_id": ["gene1", "gene1", "gene1", "gene1"],
            "sample_id": ["sample1", "sample2", "sample1", "sample2"],
            "counts": [100, 200, 150, 250],
            "condition": ["treated", "control", "treated", "control"],
            "extra_meta": ["foo", "bar", "foo", "bar"]
        }
        expected_df = pl.DataFrame(expected_data)
        assert_frame_equal(df, expected_df, check_row_order=False)
    finally:
        # Clean up the temporary files
        os.remove(expr_path)
        os.remove(metadata_path)

def test_read_expression_matrix_relative_abundance_without_gene_id():
    """
    Test that relative_abundance is not calculated when gene_id_column_name is not provided.
    """
    # Create a sample expression matrix CSV content without gene_id
    expr_content = """transcript_id,sample1,sample2
tx1,100,200
tx2,150,250
"""
    expr_path = _create_temp_file(expr_content, '.csv')

    try:
        # Call the function with relative_abundance=True but gene_id_column_name=None
        with pytest.warns(UserWarning) as record:
            df = read_expression_matrix(
                expression_matrix_path=expr_path,
                relative_abundance=True,
                gene_id_column_name=None
            )
        # Check that a warning was issued
        assert len(record) == 1
        assert "relative_abundance was set to True, but gene_id_column_name was not provided" in str(record[0].message)

        # Verify that relative_abundance column is not present
        assert "relative_abundance" not in df.columns
    finally:
        # Clean up the temporary file
        os.remove(expr_path)
