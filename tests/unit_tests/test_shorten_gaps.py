# test_shorten_gaps.py

import pytest
import polars as pl
from RNApysoforms import shorten_gaps

def test_shorten_gaps_simple_input():
    """
    Test shorten_gaps with a simple input DataFrame containing exons for a single transcript.
    """
    # Create a simple DataFrame with exons
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1"],
        "start": [100, 200, 500],
        "end": [150, 250, 600],
        "type": ["exon", "exon", "exon"],
        "strand": ["+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1"],
        "exon_number": [1, 2, 3]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that rescaled_start and rescaled_end columns are present
    assert "rescaled_start" in shortened_df.columns
    assert "rescaled_end" in shortened_df.columns

    # Check that the lengths of exons remain the same after rescaling
    original_lengths = (df["end"] - df["start"] + 1).to_list()
    rescaled_lengths = shortened_df.filter(pl.col("type") == "exon")
    rescaled_lengths = (rescaled_lengths["rescaled_end"] - rescaled_lengths["rescaled_start"] + 1).to_list()
    assert original_lengths == rescaled_lengths, "Exon lengths should remain unchanged after rescaling."

    # Check that gaps have been shortened appropriately
    # Calculate original gaps
    original_gaps = [df["start"][i+1] - df["end"][i] - 1 for i in range(len(df)-1)]
    # Calculate rescaled gaps
    rescaled_df = shortened_df.filter(pl.col("type") == "exon").sort("rescaled_start")
    rescaled_starts = rescaled_df["rescaled_start"].to_list()
    rescaled_ends = rescaled_df["rescaled_end"].to_list()
    rescaled_gaps = [rescaled_starts[i+1] - rescaled_ends[i] - 1 for i in range(len(rescaled_starts)-1)]

    # Default target_gap_width is 100
    target_gap_width = 100
    for i, gap in enumerate(rescaled_gaps):
        original_gap = original_gaps[i]
        expected_gap = min(original_gap, target_gap_width)
        assert gap == expected_gap, f"Gap {i+1} should be {expected_gap}, got {gap}."

def test_shorten_gaps_with_introns():
    """
    Test shorten_gaps with a DataFrame that already contains intron entries.
    """
    # Create a DataFrame with exons and introns
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1", "tx1", "tx1"],
        "start": [100, 151, 200, 251, 300],
        "end": [150, 199, 250, 299, 350],
        "type": ["exon", "intron", "exon", "intron", "exon"],
        "strand": ["+", "+", "+", "+", "+"],
        "seqnames": ["chr1"] * 5,
        "exon_number": [1, None, 2, None, 3]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that rescaled coordinates are present
    assert "rescaled_start" in shortened_df.columns
    assert "rescaled_end" in shortened_df.columns

    # Check that intron widths have been adjusted
    introns = shortened_df.filter(pl.col("type") == "intron")
    for width in (introns["rescaled_end"] - introns["rescaled_start"] + 1).to_list():
        assert width <= 100, "Intron widths should not exceed the target_gap_width."

def test_shorten_gaps_multiple_transcripts():
    """
    Test shorten_gaps with multiple transcripts to ensure independent rescaling.
    """
    # Create a DataFrame with multiple transcripts
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "start": [100, 500, 200, 600],
        "end": [150, 550, 250, 650],
        "type": ["exon", "exon", "exon", "exon"],
        "strand": ["+", "+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1", "chr1"],
        "exon_number": [1, 2, 1, 2]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that each transcript is rescaled independently
    for tx_id in df["transcript_id"].unique():
        tx_df = shortened_df.filter(pl.col("transcript_id") == tx_id)
        rescaled_starts = tx_df["rescaled_start"].to_list()
        rescaled_ends = tx_df["rescaled_end"].to_list()
        # Check that the exons are ordered correctly
        assert rescaled_starts == sorted(rescaled_starts), f"Exons in {tx_id} are not correctly ordered after rescaling."

def test_shorten_gaps_target_gap_width():
    """
    Test shorten_gaps with different target_gap_width values.
    """
    # Create a DataFrame with large gaps
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 1000],
        "end": [200, 1100],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 2]
    })

    # Use a custom target_gap_width
    target_gap_width = 50
    shortened_df = shorten_gaps(df, target_gap_width=target_gap_width)

    # Check that the gap has been shortened to the target_gap_width
    rescaled_df = shortened_df.filter(pl.col("type") == "exon").sort("rescaled_start")
    rescaled_starts = rescaled_df["rescaled_start"].to_list()
    rescaled_ends = rescaled_df["rescaled_end"].to_list()
    rescaled_gap = rescaled_starts[1] - rescaled_ends[0] - 1
    assert rescaled_gap == target_gap_width, f"Gap should be {target_gap_width}, got {rescaled_gap}."

def test_shorten_gaps_invalid_input_type():
    """
    Test that shorten_gaps raises a TypeError when input is not a Polars DataFrame.
    """
    # Input is a dictionary instead of a Polars DataFrame
    df = {
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"],
        "exon_number": [1]
    }

    with pytest.raises(TypeError):
        shorten_gaps(df)

def test_shorten_gaps_missing_required_columns():
    """
    Test that shorten_gaps raises a ValueError when required columns are missing.
    """
    # Missing 'exon_number' column
    df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })

    with pytest.raises(ValueError):
        shorten_gaps(df)

def test_shorten_gaps_no_gaps():
    """
    Test shorten_gaps with exons that are adjacent (no gaps) to ensure no shortening occurs.
    """
    # Exons are adjacent
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 201],
        "end": [200, 300],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 2]
    })

    shortened_df = shorten_gaps(df)

    # Check that rescaled coordinates match original lengths
    original_lengths = (df["end"] - df["start"] + 1).to_list()
    rescaled_lengths = (shortened_df["rescaled_end"] - shortened_df["rescaled_start"] + 1).to_list()
    assert original_lengths == rescaled_lengths, "Exon lengths should remain unchanged after rescaling."

    # Check that there are no gaps between rescaled exons
    rescaled_starts = shortened_df["rescaled_start"].to_list()
    rescaled_ends = shortened_df["rescaled_end"].to_list()
    rescaled_gap = rescaled_starts[1] - rescaled_ends[0] - 1
    assert rescaled_gap == 0, f"No gap expected between exons, got gap of {rescaled_gap}."


def test_shorten_gaps_different_chromosomes():
    """
    Test shorten_gaps with exons on different chromosomes to ensure it raises a ValueError.
    """
    # Exons on different chromosomes
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr2"],
        "exon_number": [1, 2]
    })

    with pytest.raises(ValueError):
        shorten_gaps(df)

def test_shorten_gaps_large_gap():
    """
    Test shorten_gaps with a very large gap to ensure proper shortening.
    """
    # Exons with a large gap between them
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 10000],
        "end": [200, 10100],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 2]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that the large gap has been shortened
    rescaled_df = shortened_df.filter(pl.col("type") == "exon").sort("rescaled_start")
    rescaled_gap = rescaled_df["rescaled_start"][1] - rescaled_df["rescaled_end"][0] - 1
    assert rescaled_gap == 100, f"Large gap should be shortened to 100, got {rescaled_gap}."

def test_shorten_gaps_zero_gap():
    """
    Test shorten_gaps with exons that are directly adjacent (zero gap).
    """
    # Exons with zero gap between them
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 201],
        "end": [200, 300],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 2]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that there is no gap after rescaling
    rescaled_df = shortened_df.filter(pl.col("type") == "exon").sort("rescaled_start")
    rescaled_gap = rescaled_df["rescaled_start"][1] - rescaled_df["rescaled_end"][0] - 1
    assert rescaled_gap == 0, f"No gap expected between exons, got {rescaled_gap}."

def test_shorten_gaps_non_standard_chromosomes():
    """
    Test shorten_gaps with non-standard chromosome names.
    """
    # Exons on non-standard chromosomes
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chrUn_gl000220", "chrUn_gl000220"],
        "exon_number": [1, 2]
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that rescaling occurred without errors
    assert len(shortened_df) == 3, "Expected 4 entries after rescaling."

def test_shorten_gaps_custom_transcript_id_column():
    """
    Test shorten_gaps with a custom transcript_id_column.
    """
    # Exons with a custom transcript ID column
    df = pl.DataFrame({
        "gene_id": ["g1", "g1"],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 2]
    })

    # Call shorten_gaps function with custom transcript_id_column
    shortened_df = shorten_gaps(df, transcript_id_column="gene_id")

    # Check that rescaling occurred correctly
    assert "rescaled_start" in shortened_df.columns
    assert "rescaled_end" in shortened_df.columns

def test_shorten_gaps_large_number_of_exons():
    """
    Test shorten_gaps with a large number of exons to assess performance and correctness.
    """
    # Create a DataFrame with 100 exons
    exon_numbers = list(range(1, 101))
    starts = list(range(100, 1000, 9))
    ends = [s + 5 for s in starts]
    df = pl.DataFrame({
        "transcript_id": ["tx1"] * 100,
        "start": starts,
        "end": ends,
        "type": ["exon"] * 100,
        "strand": ["+"] * 100,
        "seqnames": ["chr1"] * 100,
        "exon_number": exon_numbers
    })

    # Call shorten_gaps function
    shortened_df = shorten_gaps(df)

    # Check that the number of exons remains the same after rescaling
    assert len(shortened_df.filter(pl.col("type") == "exon")) == 100, "Expected 100 exons after rescaling."

    # Check that rescaled coordinates are present
    assert "rescaled_start" in shortened_df.columns
    assert "rescaled_end" in shortened_df.columns

def test_shorten_gaps_strand_specific():
    """
    Test shorten_gaps to ensure that strand information is correctly handled.
    """
    # Exons on positive and negative strands
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "exon"],
        "strand": ["+", "-"],
        "seqnames": ["chr1", "chr1"],
        "exon_number": [1, 1]
    })

    ## Transcripts from different strand should raise error
    with pytest.raises(ValueError):
        shorten_gaps(df)