# test_to_intron.py

import pytest
import polars as pl
from RNApysoforms import to_intron

def test_to_intron_simple_input():
    """
    Test to_intron with a simple input DataFrame containing exons for a single transcript.
    """
    # Create a simple DataFrame with exons
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1"],
        "start": [100, 200, 300],
        "end": [150, 250, 350],
        "type": ["exon", "exon", "exon"],
        "transcript_id": ["tx1", "tx1", "tx1"],
        "strand": ["+", "+", "+"],
        "exon_number": [1, 2, 3]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron entries have been added
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 2, "Expected 2 intron entries."

    # Check that intron positions are correct
    expected_starts = [151, 251]
    expected_ends = [199, 299]
    assert introns["start"].to_list() == expected_starts, f"Expected intron starts {expected_starts}, got {introns['start'].to_list()}."
    assert introns["end"].to_list() == expected_ends, f"Expected intron ends {expected_ends}, got {introns['end'].to_list()}."

def test_to_intron_existing_introns():
    """
    Test to_intron with a DataFrame that already contains intron entries.
    """
    # Create a DataFrame with exons and introns
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1", "chr1", "chr1"],
        "start": [100, 151, 200, 251, 300],
        "end": [150, 199, 250, 299, 350],
        "type": ["exon", "intron", "exon", "intron", "exon"],
        "transcript_id": ["tx1"] * 5,
        "strand": ["+"] * 5,
        "exon_number": [1, None, 2, None, 3]
    })

    # Call to_intron function and raise error if introns already exist
    with pytest.raises(ValueError):
        to_intron(df)

def test_to_intron_multiple_transcripts():
    """
    Test to_intron with multiple transcripts to ensure introns are calculated per transcript.
    """
    # Create a DataFrame with multiple transcripts
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1", "chr1"],
        "start": [100, 200, 150, 250],
        "end": [150, 250, 200, 300],
        "type": ["exon"] * 4,
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "strand": ["+", "+", "+", "+"],
        "exon_number": [1, 2, 1, 2]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron entries have been added for each transcript
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 2, "Expected 2 intron entries for both transcripts."

    # Check intron positions for tx1
    tx1_introns = introns.filter(pl.col("transcript_id") == "tx1")
    assert len(tx1_introns) == 1, "Expected 1 intron for tx1."
    assert tx1_introns["start"][0] == 151, "Incorrect intron start for tx1."
    assert tx1_introns["end"][0] == 199, "Incorrect intron end for tx1."

    # Check intron positions for tx2
    tx2_introns = introns.filter(pl.col("transcript_id") == "tx2")
    assert len(tx2_introns) == 1, "Expected 1 intron for tx2."
    assert tx2_introns["start"][0] == 201, "Incorrect intron start for tx2."
    assert tx2_introns["end"][0] == 249, "Incorrect intron end for tx2."

def test_to_intron_single_exon_transcript():
    """
    Test to_intron with a transcript that has only one exon (no introns should be added).
    """
    # Create a DataFrame with a single exon transcript
    df = pl.DataFrame({
        "seqnames": ["chr1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "transcript_id": ["tx1"],
        "strand": ["+"],
        "exon_number": [1]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that no intron entries have been added
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 0, "Expected no intron entries for a single exon transcript."

def test_to_intron_negative_strand():
    """
    Test to_intron with exons on the negative strand to ensure exon numbers are adjusted correctly.
    """
    # Exons on the negative strand
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [300, 100],
        "end": [350, 150],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx1"],
        "strand": ["-", "-"],
        "exon_number": [2, 1]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron entries have been added
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 1, "Expected 1 intron entry."

    # Check that exon_number for intron is adjusted correctly
    expected_exon_number = 0  # Since exon_number is decremented by 1 for negative strand
    assert introns["exon_number"][0] == expected_exon_number, f"Expected exon_number {expected_exon_number}, got {introns['exon_number'][0]}."

def test_to_intron_missing_required_columns():
    """
    Test that to_intron raises a ValueError when required columns are missing.
    """
    # Missing 'exon_number' column
    df = pl.DataFrame({
        "seqnames": ["chr1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "transcript_id": ["tx1"]
    })

    with pytest.raises(ValueError):
        to_intron(df)


def test_to_intron_auto_calculate_exon_number():
    """
    Test that to_intron automatically calculates exon_number when it's missing.
    """
    # DataFrame without 'exon_number' column
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1"],
        "start": [100, 200, 300],
        "end": [150, 250, 350],
        "type": ["exon", "exon", "exon"],
        "transcript_id": ["tx1", "tx1", "tx1"],
        "strand": ["+", "+", "+"]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that exon_number column exists in the result
    assert "exon_number" in result_df.columns, "exon_number column should be automatically calculated"
    
    # Check that exon_number values are correctly assigned
    exons = result_df.filter(pl.col("type") == "exon").sort("start")
    expected_exon_numbers = [1, 2, 3]
    assert exons["exon_number"].to_list() == expected_exon_numbers, f"Expected exon numbers {expected_exon_numbers}, got {exons['exon_number'].to_list()}."

def test_to_intron_invalid_input_type():
    """
    Test that to_intron raises a TypeError when input is not a Polars DataFrame.
    """
    # Input is a dictionary instead of a Polars DataFrame
    df = {
        "seqnames": ["chr1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "transcript_id": ["tx1"],
        "strand": ["+"],
        "exon_number": [1]
    }

    with pytest.raises(TypeError):
        to_intron(df)


def test_to_intron_exons_out_of_order():
    """
    Test to_intron with exons that are not in order to ensure sorting works correctly.
    """
    # Exons out of order
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1"],
        "start": [300, 100, 200],
        "end": [350, 150, 250],
        "type": ["exon", "exon", "exon"],
        "transcript_id": ["tx1", "tx1", "tx1"],
        "strand": ["+", "+", "+"],
        "exon_number": [3, 1, 2]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron positions are calculated correctly after sorting
    introns = result_df.filter(pl.col("type") == "intron").sort("start")
    expected_starts = [151, 251]
    expected_ends = [199, 299]
    assert introns["start"].to_list() == expected_starts, f"Expected intron starts {expected_starts}, got {introns['start'].to_list()}."
    assert introns["end"].to_list() == expected_ends, f"Expected intron ends {expected_ends}, got {introns['end'].to_list()}."

def test_to_intron_additional_features():
    """
    Test to_intron with additional genomic features like CDS to ensure they are retained.
    """
    # Exons and CDS features
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1", "chr1", "chr1"],
        "start": [100, 200, 120, 220],
        "end": [150, 250, 170, 270],
        "type": ["exon", "exon", "CDS", "CDS"],
        "transcript_id": ["tx1"] * 4,
        "strand": ["+"] * 4,
        "exon_number": [1, 2, None, None]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that CDS entries are retained
    cds_entries = result_df.filter(pl.col("type") == "CDS")
    assert len(cds_entries) == 2, "Expected 2 CDS entries."

    # Check that intron entries have been added
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 1, "Expected 1 intron entry."

def test_to_intron_non_standard_chromosomes():
    """
    Test to_intron with non-standard chromosome names.
    """
    # Exons on non-standard chromosomes
    df = pl.DataFrame({
        "seqnames": ["chrUn_gl000220", "chrUn_gl000220"],
        "start": [100, 500],
        "end": [200, 600],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx1"],
        "strand": ["+", "+"],
        "exon_number": [1, 2]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron entries have been added without errors
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 1, "Expected 1 intron entry."

def test_to_intron_custom_transcript_id_column():
    """
    Test to_intron with a custom transcript_id_column.
    """
    # Exons with a custom transcript ID column
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "exon"],
        "gene_id": ["g1", "g1"],
        "strand": ["+", "+"],
        "exon_number": [1, 2]
    })

    # Call to_intron function with custom transcript_id_column
    result_df = to_intron(df, transcript_id_column="gene_id")

    # Check that intron entries have been added
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 1, "Expected 1 intron entry with custom transcript_id_column."

def test_to_intron_large_number_of_exons():
    """
    Test to_intron with a large number of exons to assess performance and correctness.
    """
    # Create a DataFrame with 100 exons
    exon_numbers = list(range(1, 101))
    starts = list(range(100, 1000, 9))
    ends = [s + 5 for s in starts]
    df = pl.DataFrame({
        "seqnames": ["chr1"] * 100,
        "start": starts,
        "end": ends,
        "type": ["exon"] * 100,
        "transcript_id": ["tx1"] * 100,
        "strand": ["+"] * 100,
        "exon_number": exon_numbers
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that the number of introns is one less than the number of exons
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 99, "Expected 99 intron entries for 100 exons."

def test_to_intron_strand_specific():
    """
    Test to_intron to ensure that strand information is correctly handled.
    """
    # Exons on positive and negative strands
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx2"],
        "strand": ["+", "-"],
        "exon_number": [1, 1]
    })

    # Call to_intron function
    result_df = to_intron(df)

    # Check that intron entries are handled per strand
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 0, "Expected no intron entries for single exon transcripts."

def test_to_intron_no_exon_entries():
    """
    Test to_intron with a DataFrame that has no exon entries.
    """
    # DataFrame without exon entries
    df = pl.DataFrame({
        "seqnames": ["chr1"],
        "start": [100],
        "end": [200],
        "type": ["CDS"],
        "transcript_id": ["tx1"],
        "strand": ["+"],
        "exon_number": [None]
    })

    # Call to_intron function and have ValueError cause no exons
    with pytest.raises(ValueError):
        to_intron(df)


def test_to_intron_exons_with_same_start_end():
    """
    Test to_intron with exons that have the same start and end positions.
    """
    # Exons with same start and end
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [100, 100],
        "end": [200, 200],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx1"],
        "strand": ["+", "+"],
        "exon_number": [1, 2]
    })

    # Call to_intron function and have ValueError cause exons overlap
    with pytest.raises(ValueError):
        to_intron(df)
        
def test_to_intron_zero_length_introns():
    """
    Test to_intron to ensure that introns with zero or negative length are filtered out.
    """
    # Exons that result in zero-length intron
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [100, 200],
        "end": [200, 200],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx1"],
        "strand": ["+", "+"],
        "exon_number": [1, 2]
    })

    # Call to_intron function and have ValueError cause exons overlap
    with pytest.raises(ValueError):
        to_intron(df)

def test_to_intron_no_transcript_id_column():
    """
    Test to_intron when the transcript_id_column is missing.
    """
    # Missing 'transcript_id' column
    df = pl.DataFrame({
        "seqnames": ["chr1"],
        "start": [100],
        "end": [200],
        "type": ["exon"],
        "strand": ["+"],
        "exon_number": [1]
    })
    
    # Call to_intron function and have ValueError cause exons overlap
    with pytest.raises(ValueError):
        to_intron(df)

def test_to_intron_intron_positions():
    """
    Test to_intron to ensure that intron positions are correctly calculated without adjustment.
    """
    # Exons with known positions
    df = pl.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "start": [100, 300],
        "end": [200, 400],
        "type": ["exon", "exon"],
        "transcript_id": ["tx1", "tx1"],
        "strand": ["+", "+"],
        "exon_number": [1, 2]
    })

    # Expected intron positions without any adjustment
    expected_intron_start = 201
    expected_intron_end = 299

    # Call to_intron function
    result_df = to_intron(df)

    # Check intron positions
    introns = result_df.filter(pl.col("type") == "intron")
    assert len(introns) == 1, "Expected 1 intron entry."
    assert introns["start"][0] == expected_intron_start, f"Expected intron start {expected_intron_start}, got {introns['start'][0]}."
    assert introns["end"][0] == expected_intron_end, f"Expected intron end {expected_intron_end}, got {introns['end'][0]}."
