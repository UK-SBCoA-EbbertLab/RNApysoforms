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
    assert result.to_dict(as_series=False) == expected.to_dict(as_series=False)


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
    assert "The DataFrame is missing the following required columns:" in str(exc_info.value)


def test_empty_dataframe():
    empty_df = pl.DataFrame({
        "transcript_id": [],
        "start": [],
        "end": [],
        "type": [],
        "strand": []
    })

    result = calculate_exon_number(empty_df)

    print(result.head())

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
    # Annotations with multiple transcripts, exons, introns, and CDS regions
    df = pl.DataFrame({
        "transcript_id": [
            # tx1 (positive strand)
            "tx1", "tx1", "tx1", "tx1", "tx1", "tx1", "tx1",
            # tx2 (negative strand)
            "tx2", "tx2", "tx2", "tx2", "tx2", "tx2", "tx2",
            # tx3 (positive strand)
            "tx3", "tx3", "tx3", "tx3", "tx3", "tx3", "tx3"
        ],
        "start": [
            # tx1
            100, 150, 201, 301, 301, 501, 501,
            # tx2 (negative strand)
            400, 450, 601, 701, 701, 801, 801,
            # tx3
            900, 950, 1001, 1101, 1101, 1301, 1301
        ],
        "end": [
            # tx1
            200, 200, 300, 400, 400, 600, 550,
            # tx2
            500, 500, 700, 800, 800, 900, 850,
            # tx3
            1000, 1000, 1100, 1200, 1200, 1400, 1350
        ],
        "type": [
            # tx1
            "exon",    # Exon 1: 100-200
            "CDS",     # CDS 1: 150-200 (within Exon 1)
            "intron",  # Intron 1: 201-300
            "exon",    # Exon 2: 301-400
            "CDS",     # CDS 2: 301-400 (within Exon 2)
            "exon",    # Exon 3: 501-600
            "CDS",     # CDS 3: 501-550 (within Exon 3)
            # tx2
            "exon",    # Exon 1: 400-500
            "CDS",     # CDS 3: 450-500 (within Exon 1)
            "intron",  # Intron 1: 501-600
            "exon",    # Exon 2: 601-700
            "CDS",     # CDS 2: 601-700 (within Exon 2)
            "exon",    # Exon 3: 701-900
            "CDS",     # CDS 1: 801-850 (within Exon 3)
            # tx3
            "exon",    # Exon 1: 900-1000
            "CDS",     # CDS 1: 950-1000 (within Exon 1)
            "intron",  # Intron 1: 1001-1100
            "exon",    # Exon 2: 1101-1200
            "CDS",     # CDS 2: 1101-1200 (within Exon 2)
            "exon",    # Exon 3: 1301-1400
            "CDS"      # CDS 3: 1301-1350 (within Exon 3)
        ],
        "strand": [
            # tx1
            "+", "+", "+", "+", "+", "+", "+",
            # tx2
            "-", "-", "-", "-", "-", "-", "-",
            # tx3
            "+", "+", "+", "+", "+", "+", "+"
        ]
    })

    # Expected exon numbers for each row
    expected_exon_numbers = [
        # tx1
        1, 1, 1, 2, 2, 3, 3,
        # tx2 (negative strand)
        3, 3, 3, 2, 2, 1, 1,
        # tx3
        1, 1, 1, 2, 2, 3, 3
    ]

    # Call the function to calculate exon numbers
    result = calculate_exon_number(df)

    # Assert that the calculated exon numbers match the expected values
    assert result["exon_number"].to_list() == expected_exon_numbers, \
        f"Expected {expected_exon_numbers}, but got {result['exon_number'].to_list()}"