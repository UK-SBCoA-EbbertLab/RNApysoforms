import os
import polars as pl
import plotly.graph_objects as go
import RNApysoforms as RNApy
import pytest

def test_wonky_gtf_visualization():
    """
    Test end-to-end functionality for processing and visualizing a custom GTF file.
    
    This test demonstrates:
    1. Reading and parsing a custom GTF file
    2. Extracting attributes from the GTF
    3. Calculating exon numbers
    4. Shortening intron gaps for visualization
    5. Creating and displaying transcript structure traces
    """
    # Get the directory of the current file
    current_dir = os.path.dirname(__file__)
    # Construct the path to the test_data directory
    test_data_dir = os.path.abspath(os.path.join(current_dir, '../test_data/'))
    
    # Define the path to the custom GTF file
    gtf_path = os.path.join(test_data_dir, "MYB.annot.gtf")
    
    # Define the column names for the GTF file (standard GTF format)
    column_names = [
        "seqnames",    # Chromosome or sequence name
        "source",      # Annotation source
        "type",        # Feature type (e.g., exon, CDS)
        "start",       # Start position of the feature
        "end",         # End position of the feature
        "score",       # Score value (usually '.')
        "strand",      # Strand information ('+' or '-')
        "phase",       # Reading frame phase
        "attributes"   # Additional attributes in key-value pairs
    ]

    # Define data types for GTF columns
    dtypes = {
        "seqnames": pl.Utf8,
        "source": pl.Utf8,
        "type": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "phase": pl.Utf8,
        "attributes": pl.Utf8
    }

    # Read the GTF file using Polars
    gtf = pl.read_csv(
        gtf_path,
        separator="\t",
        has_header=False,
        comment_prefix="#",       # Skip comment lines starting with '#'
        new_columns=column_names, # Assign column names since GTF files have no header
        schema_overrides=dtypes   # Specify data types for each column
    )

    # Remove unnecessary columns (phase and score)
    gtf = gtf.select([col for col in column_names if col not in ["phase", "score"]])

    # Extract attributes from the GTF file and create new columns
    gtf = gtf.with_columns([
        pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),
        pl.col("attributes").str.extract(r'transcript_status "([^"]+)"', 1).alias("transcript_status"),
        pl.col("attributes").str.extract(r'transcript_type "([^"]+)"', 1).alias("transcript_type"),
        pl.col("attributes").str.extract(r'transcript_category "([^"]+)"', 1).alias("transcript_category"),
        pl.col("attributes").str.extract(r'CDS_id "([^"]+)"', 1).alias("CDS_id"),
        pl.col("attributes").str.extract(r'CDS_status "([^"]+)"', 1).alias("CDS_status")
    ])
    
    # Filter to specific transcripts of interest
    transcripts_to_keep = [
        "MYB:::ENST00000420123.6__NA__remap-9569__NA", 
        "MYB:::NA__PB.7027.4__NA__NA",
        "MYB:::NA__PB.7027.3__NA__NA", 
        "MYB:::ENST00000367814.8__PB.7027.1__NA__NA"
    ]
    filtered_gtf = gtf.filter(pl.col("transcript_id").is_in(transcripts_to_keep))

    # Shorten intron gaps for better visualization
    shortened_gtf = RNApy.shorten_gaps(filtered_gtf)

    # Create visualization traces
    traces = RNApy.make_traces(
        annotation=shortened_gtf, 
        x_start="rescaled_start", 
        x_end="rescaled_end",
        y='transcript_id', 
        annotation_hue="CDS_status"
    )
    
    # Generate the plot
    fig = RNApy.make_plot(
        traces=traces, 
        subplot_titles=["Transcript Structure"], 
        width=900, 
        height=500
    )

    # Verify that the figure was created correctly
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0
    assert "Transcript Structure" in fig.layout.annotations[0].text
