# tests/test_gene_filtering.py

import pytest
import polars as pl
import warnings
from RNApysoforms import gene_filtering

def test_gene_filtering_only_annotation():
    """
    Test filtering with only the annotation DataFrame provided.
    """
    # Create a sample annotation DataFrame
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1", "TP53"],
        "transcript_id": ["tx1", "tx2", "tx3"],
        "counts": [100, 200, 150]
    })
    target_gene = "BRCA1"

    # Call the function
    filtered_annotation = gene_filtering(target_gene, annotation_df)

    # Assert that the filtered_annotation contains only entries for BRCA1
    assert filtered_annotation["gene_name"].unique().to_list() == ["BRCA1"]
    # Assert that the number of entries is correct
    assert filtered_annotation.shape[0] == 2
    # Check that the transcripts are tx1 and tx2
    assert set(filtered_annotation["transcript_id"].to_list()) == {"tx1", "tx2"}

def test_gene_filtering_with_expression_matrix():
    """
    Test filtering with both annotation and expression matrix provided.
    """
    # Create sample annotation DataFrame
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1", "TP53", "BRCA1"],
        "transcript_id": ["tx1", "tx2", "tx3", "tx4"],
        "other_info": [1, 2, 3, 4]
    })

    # Create sample expression matrix DataFrame
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2", "tx4", "tx5"],
        "counts": [100, 200, 300, 400],
        "sample_id": ["sample1", "sample2", "sample3", "sample4"]
    })

    target_gene = "BRCA1"

    # Call the function
    filtered_annotation, filtered_expression_matrix = gene_filtering(
        target_gene,
        annotation_df,
        expression_matrix=expression_matrix_df
    )

    # Assert that the filtered_annotation contains only entries for BRCA1
    assert filtered_annotation["gene_name"].unique().to_list() == ["BRCA1"]
    # Assert that the number of entries is correct
    assert filtered_annotation.shape[0] == 3
    # Assert that the filtered_expression_matrix contains the correct transcripts
    assert set(filtered_expression_matrix["transcript_id"].unique().to_list()) == {"tx1", "tx2", "tx4"}

def test_gene_filtering_no_matching_gene():
    """
    Test that ValueError is raised when the target gene is not found in the annotation.
    """
    # Create a sample annotation DataFrame
    annotation_df = pl.DataFrame({
        "gene_name": ["TP53", "EGFR"],
        "transcript_id": ["tx1", "tx2"]
    })
    target_gene = "BRCA1"

    # Call the function and expect a ValueError
    with pytest.raises(ValueError) as excinfo:
        gene_filtering(target_gene, annotation_df)
    assert f"No annotation found for gene: {target_gene}" in str(excinfo.value)

def test_gene_filtering_invalid_annotation_type():
    """
    Test that TypeError is raised when the annotation is not a Polars DataFrame.
    """
    annotation_df = {"gene_name": ["BRCA1"], "transcript_id": ["tx1"]}
    target_gene = "BRCA1"

    with pytest.raises(TypeError) as excinfo:
        gene_filtering(target_gene, annotation_df)
    assert "Expected 'annotation' to be of type pl.DataFrame" in str(excinfo.value)

def test_gene_filtering_missing_required_columns():
    """
    Test that ValueError is raised when required columns are missing in the annotation.
    """
    # Create a sample annotation DataFrame missing the 'gene_name' column
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "counts": [100, 200]
    })
    target_gene = "BRCA1"

    with pytest.raises(ValueError) as excinfo:
        gene_filtering(target_gene, annotation_df)
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_gene_filtering_expression_matrix_empty_after_filtering():
    """
    Test that ValueError is raised when the expression matrix is empty after filtering.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1"],
        "transcript_id": ["tx1"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx2"],
        "counts": [100]
    })
    target_gene = "BRCA1"

    with pytest.raises(ValueError) as excinfo:
        gene_filtering(target_gene, annotation_df, expression_matrix=expression_matrix_df)
    assert "Expression matrix is empty after filtering" in str(excinfo.value)

def test_gene_filtering_warning_missing_transcripts_in_expression():
    """
    Test that a warning is issued when transcripts are missing in the expression matrix.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1"],
        "transcript_id": ["tx1", "tx2"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "counts": [100]
    })
    target_gene = "BRCA1"

    with pytest.warns(UserWarning) as record:
        filtered_annotation, filtered_expression_matrix = gene_filtering(
            target_gene,
            annotation_df,
            expression_matrix=expression_matrix_df
        )

    assert len(record) == 1
    assert "transcript(s) are present in the annotation but missing in the expression matrix" in str(record[0].message)

def test_gene_filtering_order_by_expression():
    """
    Test ordering of transcripts by expression levels.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1", "BRCA1"],
        "transcript_id": ["tx1", "tx2", "tx3"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2", "tx3"],
        "counts": [300, 100, 200]
    })
    target_gene = "BRCA1"

    filtered_annotation, filtered_expression_matrix = gene_filtering(
        target_gene,
        annotation_df,
        expression_matrix=expression_matrix_df,
        order_by_expression=True
    )

    # Assert that transcripts are ordered by expression counts descending
    expected_order = ["tx2", "tx3", "tx1"]
    assert filtered_annotation["transcript_id"].to_list() == expected_order
    assert filtered_expression_matrix["transcript_id"].unique(maintain_order=True).to_list() == expected_order

def test_gene_filtering_keep_top_expressed_transcripts():
    """
    Test keeping only the top N expressed transcripts.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1", "BRCA1", "BRCA1"],
        "transcript_id": ["tx1", "tx2", "tx3", "tx4"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2", "tx3", "tx4"],
        "counts": [400, 300, 200, 100]
    })
    target_gene = "BRCA1"

    filtered_annotation, filtered_expression_matrix = gene_filtering(
        target_gene,
        annotation_df,
        expression_matrix=expression_matrix_df,
        keep_top_expressed_transcripts=2
    )

    # Assert that only top 2 transcripts are kept
    expected_transcripts = ["tx1", "tx2"]
    assert set(filtered_annotation["transcript_id"].to_list()) == set(expected_transcripts)
    assert set(filtered_expression_matrix["transcript_id"].unique().to_list()) == set(expected_transcripts)

def test_gene_filtering_invalid_keep_top_expressed_transcripts():
    """
    Test that ValueError is raised when 'keep_top_expressed_transcripts' is invalid.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1"],
        "transcript_id": ["tx1"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "counts": [100]
    })
    target_gene = "BRCA1"

    with pytest.raises(ValueError) as excinfo:
        gene_filtering(
            target_gene,
            annotation_df,
            expression_matrix=expression_matrix_df,
            keep_top_expressed_transcripts=0
        )
    assert "'keep_top_expressed_transcripts' must be 'all' or a positive integer" in str(excinfo.value)

def test_gene_filtering_keep_top_transcripts_exceeds_available():
    """
    Test that a warning is issued when 'keep_top_expressed_transcripts' exceeds available transcripts.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1"],
        "transcript_id": ["tx1", "tx2"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "counts": [100, 200]
    })
    target_gene = "BRCA1"

    with pytest.warns(UserWarning) as record:
        filtered_annotation, filtered_expression_matrix = gene_filtering(
            target_gene,
            annotation_df,
            expression_matrix=expression_matrix_df,
            keep_top_expressed_transcripts=5
        )

    assert len(record) == 1
    assert "exceeds the total number of transcripts" in str(record[0].message)
    # Assert that all transcripts are kept
    assert filtered_annotation.shape[0] == 2
    assert filtered_expression_matrix.shape[0] == 2

def test_gene_filtering_order_by_expression_false():
    """
    Test that transcripts are not ordered when 'order_by_expression' is False.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1", "BRCA1", "BRCA1"],
        "transcript_id": ["tx3", "tx1", "tx2"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2", "tx3"],
        "counts": [100, 200, 300]
    })
    target_gene = "BRCA1"

    filtered_annotation, filtered_expression_matrix = gene_filtering(
        target_gene,
        annotation_df,
        expression_matrix=expression_matrix_df,
        order_by_expression=False
    )

    # Assert that the order remains as in the annotation DataFrame
    expected_order = ["tx3", "tx1", "tx2"]
    assert filtered_annotation["transcript_id"].to_list() == expected_order

def test_gene_filtering_transcripts_in_expression_not_in_annotation():
    """
    Test that transcripts present in expression matrix but not in annotation are ignored without warning.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1"],
        "transcript_id": ["tx1"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "counts": [100, 200]
    })
    target_gene = "BRCA1"

    # Capture warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        filtered_annotation, filtered_expression_matrix = gene_filtering(
            target_gene,
            annotation_df,
            expression_matrix=expression_matrix_df
        )
        # Assert that no warnings were raised
        assert len(w) == 0

    # Assert that only tx1 is in the filtered_expression_matrix
    assert filtered_expression_matrix["transcript_id"].unique().to_list() == ["tx1"]

def test_gene_filtering_custom_column_names():
    """
    Test using custom column names for transcript_id and gene_name.
    """
    # Create sample annotation and expression matrix DataFrames with custom column names
    annotation_df = pl.DataFrame({
        "gene_id": ["BRCA1", "BRCA1", "TP53"],
        "tx_id": ["tx1", "tx2", "tx3"]
    })
    expression_matrix_df = pl.DataFrame({
        "tx_id": ["tx1", "tx2"],
        "expression": [100, 200]
    })
    target_gene = "BRCA1"

    filtered_annotation, filtered_expression_matrix = gene_filtering(
        target_gene,
        annotation_df,
        expression_matrix=expression_matrix_df,
        transcript_id_column="tx_id",
        gene_id_column="gene_id",
        order_by_expression_column="expression"
    )

    # Assert that the filtered_annotation contains only BRCA1 entries
    assert filtered_annotation["gene_id"].unique().to_list() == ["BRCA1"]
    # Assert that the filtered_expression_matrix contains the correct transcripts
    assert set(filtered_expression_matrix["tx_id"].unique().to_list()) == {"tx1", "tx2"}

def test_gene_filtering_no_expression_column_in_expression_matrix():
    """
    Test that ValueError is raised when the order_by_expression_column is missing in the expression matrix.
    """
    # Create sample annotation and expression matrix DataFrames
    annotation_df = pl.DataFrame({
        "gene_name": ["BRCA1"],
        "transcript_id": ["tx1"]
    })
    expression_matrix_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "other_counts": [100]
    })
    target_gene = "BRCA1"

    with pytest.raises(ValueError) as excinfo:
        gene_filtering(
            target_gene,
            annotation_df,
            expression_matrix=expression_matrix_df,
            order_by_expression_column="counts"
        )
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)
