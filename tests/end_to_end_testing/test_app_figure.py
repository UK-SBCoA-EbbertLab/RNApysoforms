import pytest
import os
import RNApysoforms as RNApy
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

def test_plot_APP_gene():
    """
    Test end-to-end functionality for visualizing the APP gene with expression data.
    
    This test demonstrates:
    1. Reading and parsing an Ensembl GTF file
    2. Loading and normalizing expression data
    3. Filtering annotation and expression data for a specific gene
    4. Shortening intron gaps for better visualization
    5. Creating transcript structure traces with expression data
    6. Generating a multi-panel figure with transcript structure and expression metrics
    
    The test verifies that the final figure is correctly created with the expected components.
    """
    # Get the directory of the current file
    current_dir = os.path.dirname(__file__)
    # Construct the path to the test_data directory
    test_data_dir = os.path.abspath(os.path.join(current_dir, '../test_data/'))

    # Define file paths
    gtf_path = os.path.join(test_data_dir, "Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")
    expression_matrix_path = os.path.join(test_data_dir, "counts_matrix_chr21_and_Y.tsv")
    metadata_path = os.path.join(test_data_dir, "sample_metadata.tsv")

    # Read GTF file
    annotation = RNApy.read_ensembl_gtf(gtf_path)

    # Read expression matrix
    counts = RNApy.read_expression_matrix(
        expression_matrix_path=expression_matrix_path,
        metadata_path=metadata_path,
        cpm_normalization=True,  # Calculate Counts Per Million
        relative_abundance=True  # Calculate relative abundance within genes
    )

    # Define a mapping from transcript_biotype to colors
    biotype_colors = {
        'protein_coding': '#F8766D',
        'retained_intron': '#00BFC4',
        'protein_coding_CDS_not_defined': 'green'
    }

    # Define gene name to filter
    gene_name = "APP"

    # Filter gene name in annotation and counts matrix
    annotation, counts = RNApy.gene_filtering(
        annotation=annotation,
        expression_matrix=counts,
        target_gene=gene_name,
        transcript_id_column="transcript_id",
        order_by_expression=True,  # Order transcripts by expression level
        keep_top_expressed_transcripts="all",  # Keep all transcripts
        order_by_expression_column="counts"  # Use raw counts for ordering
    )

    # Shorten gaps for better visualization
    rescaled_annotation = RNApy.shorten_gaps(
        annotation=annotation,
        transcript_id_column="transcript_id"
    )

    # Create traces for visualization
    traces = RNApy.make_traces(
        annotation=rescaled_annotation,
        expression_matrix=counts,
        order_transcripts_by_expression_matrix=True,
        y='transcript_id',
        expression_columns=["counts", "relative_abundance", "CPM"],  # Show multiple expression metrics
        x_start="rescaled_start",
        x_end="rescaled_end",
        annotation_fill_color="blue",
        expression_fill_color="grey",
        exon_height=0.3,
        cds_height=0.5,
        arrow_size=9,
        annotation_hue="transcript_biotype",  # Color transcripts by biotype
        expression_hue="AD status",  # Color expression by Alzheimer's Disease status
        expression_plot_style="boxplot",
        spanmode="hard"
    )

    # Generate the plot with multiple panels
    fig = RNApy.make_plot(
        traces=traces,
        subplot_titles=["Transcript Structure", "Counts", "Relative Abundance", "CPM"],
        showlegend=True,
        boxgap=0.2,
        boxgroupgap=0,
        vert_grid_transcript_structure_plot=False
    )

    # Verify that the figure was created correctly
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0
    assert "Transcript Structure" in fig.layout.annotations[0].text

