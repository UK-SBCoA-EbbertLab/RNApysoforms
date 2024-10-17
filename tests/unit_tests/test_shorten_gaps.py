# test_shorten_gaps.py

import polars as pl
import pytest
from RNApysoforms import shorten_gaps  # Replace with the actual module name where shorten_gaps is defined

def test_shorten_gaps_basic():
    """
    Test the basic functionality with a simple DataFrame containing exons.
    Verify that introns are generated and included in the output.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1"],
        "exon_number": [1, 2, 3],
        "start": [100, 300, 600],
        "end": [200, 400, 700],
        "type": ["exon", "exon", "exon"],
        "strand": ["+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that introns are included in the output
    assert "intron" in shortened_df["type"].unique().to_list()

    # Separate exons and introns for verification
    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")
    introns_df = shortened_df.filter(pl.col("type") == "intron").sort("exon_number")

    # Expected rescaled coordinates for exons
    expected_exon_rescaled_starts = [2, 202, 502]
    expected_exon_rescaled_ends = [102, 302, 602]

    # Expected rescaled coordinates for introns
    expected_intron_rescaled_starts = [102, 302]
    expected_intron_rescaled_ends = [202, 502]

    # Verify exon rescaled coordinates
    assert exons_df["rescaled_start"].to_list() == expected_exon_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_exon_rescaled_ends

    # Verify intron rescaled coordinates
    assert introns_df["rescaled_start"].to_list() == expected_intron_rescaled_starts
    assert introns_df["rescaled_end"].to_list() == expected_intron_rescaled_ends

def test_shorten_gaps_with_existing_introns():
    """
    Test the function when intron entries are already present in the DataFrame.
    Verify that introns are correctly processed and rescaled.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1", "tx1", "tx1"],
        "exon_number": [1, None, 2, None, 3],
        "start": [100, 201, 300, 401, 600],
        "end": [200, 300, 400, 500, 700],
        "type": ["exon", "intron", "exon", "intron", "exon"],
        "strand": ["+", "+", "+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1", "chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that introns are included in the output
    assert "intron" in shortened_df["type"].unique().to_list()

    # Separate exons and introns
    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")
    introns_df = shortened_df.filter(pl.col("type") == "intron").sort("rescaled_start")

    # Expected rescaled coordinates for exons
    expected_exon_rescaled_starts = [2, 201, 400]
    expected_exon_rescaled_ends = [102, 301, 500]

    # Expected rescaled coordinates for introns
    expected_intron_rescaled_starts = [102, 301]
    expected_intron_rescaled_ends = [201, 400]

    # Verify exon rescaled coordinates
    assert exons_df["rescaled_start"].to_list() == expected_exon_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_exon_rescaled_ends

    # Verify intron rescaled coordinates
    assert introns_df["rescaled_start"].to_list() == expected_intron_rescaled_starts
    assert introns_df["rescaled_end"].to_list() == expected_intron_rescaled_ends

def test_shorten_gaps_no_introns():
    """
    Test the function when intron entries are not present and need to be generated.
    Verify that introns are generated and correctly processed.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that introns are included in the output
    assert "intron" in shortened_df["type"].unique().to_list()

    # Separate exons and introns
    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")
    introns_df = shortened_df.filter(pl.col("type") == "intron")

    # Expected rescaled coordinates for exons
    expected_exon_rescaled_starts = [2, 402]
    expected_exon_rescaled_ends = [102, 502]

    # Expected rescaled coordinates for introns
    expected_intron_rescaled_starts = [102]
    expected_intron_rescaled_ends = [402]

    # Verify exon rescaled coordinates
    assert exons_df["rescaled_start"].to_list() == expected_exon_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_exon_rescaled_ends

    # Verify intron rescaled coordinates
    assert introns_df["rescaled_start"].to_list() == expected_intron_rescaled_starts
    assert introns_df["rescaled_end"].to_list() == expected_intron_rescaled_ends

def test_shorten_gaps_multiple_transcripts():
    """
    Test the function with multiple transcripts.
    Verify that each transcript is processed independently and introns are correctly handled.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "exon_number": [1, 2, 1, 2],
        "start": [100, 500, 200, 600],
        "end": [200, 600, 300, 700],
        "type": ["exon", "exon", "exon", "exon"],
        "strand": ["+", "+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that introns are included
    assert "intron" in shortened_df["type"].unique().to_list()

    # Separate the transcripts
    tx1_df = shortened_df.filter(pl.col("transcript_id") == "tx1")
    tx2_df = shortened_df.filter(pl.col("transcript_id") == "tx2")

    # Expected rescaled coordinates for tx1 exons
    tx1_exons = tx1_df.filter(pl.col("type") == "exon").sort("exon_number")
    expected_tx1_exon_rescaled_starts = [2, 253]
    expected_tx1_exon_rescaled_ends = [102, 353]

    # Expected rescaled coordinates for tx1 introns
    tx1_introns = tx1_df.filter(pl.col("type") == "intron")
    expected_tx1_intron_rescaled_starts = [102]
    expected_tx1_intron_rescaled_ends = [253]

    # Verify tx1 exons
    assert tx1_exons["rescaled_start"].to_list() == expected_tx1_exon_rescaled_starts
    assert tx1_exons["rescaled_end"].to_list() == expected_tx1_exon_rescaled_ends

    # Verify tx1 introns
    assert tx1_introns["rescaled_start"].to_list() == expected_tx1_intron_rescaled_starts
    assert tx1_introns["rescaled_end"].to_list() == expected_tx1_intron_rescaled_ends

    # Expected rescaled coordinates for tx2 exons
    tx2_exons = tx2_df.filter(pl.col("type") == "exon").sort("exon_number")
    expected_tx2_exon_rescaled_starts = [102, 353]
    expected_tx2_exon_rescaled_ends = [202, 453]

    # Expected rescaled coordinates for tx2 introns
    tx2_introns = tx2_df.filter(pl.col("type") == "intron")
    expected_tx2_intron_rescaled_starts = [202]
    expected_tx2_intron_rescaled_ends = [353]

    # Verify tx2 exons
    assert tx2_exons["rescaled_start"].to_list() == expected_tx2_exon_rescaled_starts
    assert tx2_exons["rescaled_end"].to_list() == expected_tx2_exon_rescaled_ends

    # Verify tx2 introns
    assert tx2_introns["rescaled_start"].to_list() == expected_tx2_intron_rescaled_starts
    assert tx2_introns["rescaled_end"].to_list() == expected_tx2_intron_rescaled_ends

def test_shorten_gaps_different_strands():
    """
    Test the function with exons on different strands.
    Verify that exons on different strands are processed separately.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "exon_number": [1, 2, 1, 2],
        "start": [100, 500, 200, 600],
        "end": [200, 600, 300, 700],
        "type": ["exon", "exon", "exon", "exon"],
        "strand": ["+", "+", "-", "-"],
        "seqnames": ["chr1", "chr1", "chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that introns are included
    assert "intron" in shortened_df["type"].unique().to_list()

    # Verify that tx1 and tx2 are processed separately
    tx1_df = shortened_df.filter(pl.col("transcript_id") == "tx1")
    tx2_df = shortened_df.filter(pl.col("transcript_id") == "tx2")

    # Expected rescaled coordinates for tx1 exons
    tx1_exons = tx1_df.filter(pl.col("type") == "exon").sort("exon_number")
    expected_tx1_exon_rescaled_starts = [2, 202]
    expected_tx1_exon_rescaled_ends = [102, 302]

    # Expected rescaled coordinates for tx2 exons
    tx2_exons = tx2_df.filter(pl.col("type") == "exon").sort("exon_number")
    expected_tx2_exon_rescaled_starts = [2, 202]
    expected_tx2_exon_rescaled_ends = [102, 302]

    # Verify tx1 exons
    assert tx1_exons["rescaled_start"].to_list() == expected_tx1_exon_rescaled_starts
    assert tx1_exons["rescaled_end"].to_list() == expected_tx1_exon_rescaled_ends

    # Verify tx2 exons
    assert tx2_exons["rescaled_start"].to_list() == expected_tx2_exon_rescaled_starts
    assert tx2_exons["rescaled_end"].to_list() == expected_tx2_exon_rescaled_ends

def test_shorten_gaps_missing_required_columns():
    """
    Test that the function raises a ValueError when required columns are missing.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        # Missing 'type', 'strand', 'seqnames', 'exon_number'
    })

    with pytest.raises(ValueError) as exc_info:
        shorten_gaps(df, transcript_id_column="transcript_id")

    assert "The DataFrame is missing the following required columns" in str(exc_info.value)

def test_shorten_gaps_incorrect_input_type():
    """
    Test that the function raises a TypeError when the input is not a Polars DataFrame.
    """
    df = {
        "transcript_id": ["tx1"],
        "exon_number": [1],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    }

    with pytest.raises(TypeError) as exc_info:
        shorten_gaps(df, transcript_id_column="transcript_id")

    assert "Expected 'annotation' to be of type pl.DataFrame" in str(exc_info.value)

def test_shorten_gaps_empty_dataframe():
    """
    Test the function with an empty DataFrame.
    Verify that the function returns an empty DataFrame with expected columns.
    """
    df = pl.DataFrame({
        "transcript_id": [],
        "exon_number": [],
        "start": [],
        "end": [],
        "type": [],
        "strand": [],
        "seqnames": []
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id")

    # The result should be an empty DataFrame with expected columns
    expected_columns = df.columns + ["rescaled_start", "rescaled_end"]
    assert shortened_df.shape == (0, len(expected_columns))
    assert shortened_df.columns == expected_columns

def test_shorten_gaps_negative_gap_width():
    """
    Test that the function raises an error when target_gap_width is negative.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    with pytest.raises(ValueError) as exc_info:
        shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=-50)

    assert "target_gap_width must be a positive integer" in str(exc_info.value)

def test_shorten_gaps_zero_gap_width():
    """
    Test the function with target_gap_width set to zero.
    Verify that introns are minimized to zero length.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=0)

    # Expected rescaled coordinates for exons
    expected_exon_rescaled_starts = [2, 102]
    expected_exon_rescaled_ends = [102, 202]

    # Introns minimized to zero length
    introns_df = shortened_df.filter(pl.col("type") == "intron")
    intron_lengths = (introns_df["rescaled_end"] - introns_df["rescaled_start"]).to_list()
    assert all(length == 0 for length in intron_lengths)

def test_shorten_gaps_with_cds_regions():
    """
    Test the function with CDS regions present.
    Verify that CDS regions are correctly rescaled along with exons and introns.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1"],
        "exon_number": [1, 2, None],
        "start": [100, 300, 150],
        "end": [200, 400, 250],
        "type": ["exon", "exon", "CDS"],
        "strand": ["+", "+", "+"],
        "seqnames": ["chr1", "chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Verify that CDS regions are included and rescaled
    cds_df = shortened_df.filter(pl.col("type") == "CDS")
    assert not cds_df.is_empty()

    # Expected rescaled coordinates for CDS
    expected_cds_rescaled_start = [52]
    expected_cds_rescaled_end = [152]

    assert cds_df["rescaled_start"].to_list() == expected_cds_rescaled_start
    assert cds_df["rescaled_end"].to_list() == expected_cds_rescaled_end

def test_shorten_gaps_large_target_gap_width():
    """
    Test the function with a large target_gap_width that exceeds actual gaps.
    Verify that the gaps are not altered in this case.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 150],
        "end": [110, 160],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=1000)

    # Expected rescaled coordinates (should remain proportional)
    expected_exon_rescaled_starts = [2, 52]
    expected_exon_rescaled_ends = [12, 62]

    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")

    assert exons_df["rescaled_start"].to_list() == expected_exon_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_exon_rescaled_ends

def test_shorten_gaps_transcript_id_missing():
    """
    Test the function when transcript_id_column is None.
    Verify that the function processes all exons together.
    """
    df = pl.DataFrame({
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column=None, target_gap_width=50)

    # Expected rescaled coordinates
    expected_rescaled_starts = [2, 202]
    expected_rescaled_ends = [102, 302]

    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")

    assert exons_df["rescaled_start"].to_list() == expected_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_rescaled_ends

def test_shorten_gaps_non_standard_chromosomes():
    """
    Test the function with non-standard chromosome names.
    Verify that it handles them correctly.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chrUn", "chrUn"]
    })

    shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)

    # Expected rescaled coordinates
    expected_rescaled_starts = [2, 202]
    expected_rescaled_ends = [102, 302]

    exons_df = shortened_df.filter(pl.col("type") == "exon").sort("exon_number")

    assert exons_df["rescaled_start"].to_list() == expected_rescaled_starts
    assert exons_df["rescaled_end"].to_list() == expected_rescaled_ends

def test_shorten_gaps_chromosome_mismatch():
    """
    Test that the function raises a ValueError when exons are from different chromosomes.
    """
    df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "exon_number": [1, 2],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr2"]
    })

    with pytest.raises(ValueError) as exc_info:
        shorten_gaps(df, transcript_id_column="transcript_id")

    assert "Exons must be from a single chromosome and strand" in str(exc_info.value)
