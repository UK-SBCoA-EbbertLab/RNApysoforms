import polars as pl
import pytest
from RNApysoforms import calculate_exon_number


def test_basic_functionality():
    # Create a simple annotation DataFrame
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1"],
        "start": [100, 200, 300],
        "end": [150, 250, 350],
        "type": ["exon", "exon", "exon"],
        "strand": ["+", "+", "+"]
    })

    expected = df.with_columns([
        pl.Series("exon_number", [1, 2, 3])
    ])

    result = calculate_exon_number(df)

    # Assert that the result matches the expected output
    assert result.to_dict(False) == expected.to_dict(False)


def test_strand_direction():
    # Create annotations for both positive and negative strands
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "start": [100, 200, 300, 400],
        "end": [150, 250, 350, 450],
        "type": ["exon", "exon", "exon", "exon"],
        "strand": ["+", "+", "-", "-"]
    })

    expected_exon_numbers = [1, 2, 2, 1]  # Positive strand increases, negative strand decreases

    result = calculate_exon_number(df)

    assert result.sort("transcript_id")["exon_number"].to_list() == expected_exon_numbers


def test_missing_required_columns():
    df_missing_columns = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        # 'end' column is missing
        "type": ["exon"],
        "strand": ["+"]
    })

    with pytest.raises(ValueError) as exc_info:
        calculate_exon_number(df_missing_columns)
    assert "The following required columns are missing" in str(exc_info.value)


def test_empty_dataframe():
    empty_df = pl.DataFrame({
        "transcript_id": [],
        "start": [],
        "end": [],
        "type": [],
        "strand": []
    })

    result = calculate_exon_number(empty_df)

    assert result.is_empty()


def test_non_polars_dataframe():
    import pandas as pd

    df_pandas = pd.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"]
    })

    with pytest.raises(TypeError) as exc_info:
        calculate_exon_number(df_pandas)
    assert "Expected 'annotation' to be of type pl.DataFrame" in str(exc_info.value)


def test_complex_transcript_structure():
    # Annotations with exons, introns, and CDS regions
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1", "tx1", "tx1"],
        "start": [100, 151, 201, 251, 301],
        "end": [150, 200, 250, 300, 350],
        "type": ["exon", "intron", "CDS", "exon", "CDS"],
        "strand": ["+", "+", "+", "+", "+"]
    })

    expected_exon_numbers = [1, 1, 2, 2, 3]

    result = calculate_exon_number(df)

    assert result["exon_number"].to_list() == expected_exon_numbers
