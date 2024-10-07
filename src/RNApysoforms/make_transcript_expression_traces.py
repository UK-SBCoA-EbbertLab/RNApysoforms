import plotly.graph_objects as go
import polars as pl
from typing import List, Optional, Dict

def make_transcript_expression_traces(
    long_expression_matrix: pl.DataFrame,
    y: str = "transcript_id",
    x: str = "counts",
    hue: Optional[str] = None,  # New hue parameter
    line_width: float = 0.5,
    sample_id_column: str = "sample_id",
    color_map: dict = None,  # Optional color map for hues
) -> List[go.Box]:
    """
    Create a list of Plotly Box traces for each unique transcript, colored and dodged by a hue variable.

    Parameters:
    - long_expression_matrix (pl.DataFrame): DataFrame containing transcript expression data.
    - y (str): Column name for transcript IDs (categorical axis).
    - x (str): Column name for counts (numerical axis).
    - hue (str, optional): Column name for the hue (grouping) variable.
    - line_width (float): Width of the box plot lines.
    - sample_id_column (str): Column name for sample IDs.
    - color_map (dict, optional): Dictionary mapping hue values to colors.

    Returns:
    - List[go.Box]: List of Plotly Box traces.
    """
    
    # Ensure the necessary columns exist
    required_columns = [y, x, sample_id_column]
    if hue:
        required_columns.append(hue)
    missing_columns = [col for col in required_columns if col not in long_expression_matrix.columns]
    if missing_columns:
        raise ValueError(f"Columns {missing_columns} not found in the DataFrame.")
    
    # If hue is provided, handle coloring and dodging
    if hue:
        unique_hues = long_expression_matrix[hue].unique().to_list()
        
        # Set up default color map if none provided
        if color_map is None:
            # Generate distinct colors for each hue
            default_colors = [
                'rgba(31, 119, 180, 0.6)',
                'rgba(255, 127, 14, 0.6)',
                'rgba(44, 160, 44, 0.6)',
                'rgba(214, 39, 40, 0.6)',
                'rgba(148, 103, 189, 0.6)',
                'rgba(140, 86, 75, 0.6)',
                'rgba(227, 119, 194, 0.6)',
                'rgba(127, 127, 127, 0.6)',
                'rgba(188, 189, 34, 0.6)',
                'rgba(23, 190, 207, 0.6)'
            ]
            color_map = {hue_val: default_colors[i % len(default_colors)] for i, hue_val in enumerate(unique_hues)}
        
        # Initialize list of Box traces
        traces_list = []
        
        # Iterate over each unique hue to create Box traces
        for hue_val in unique_hues:
            # Filter data for the current hue
            hue_filtered_df = long_expression_matrix.filter(pl.col(hue) == hue_val)
            
            # Create a Box trace for the current hue
            box_trace = go.Box(
                y=hue_filtered_df[y].to_list(),  # Categorical axis
                x=hue_filtered_df[x].to_list(),  # Numerical data
                text=hue_filtered_df[sample_id_column].to_list(),  # Data point labels
                name=str(hue_val),  # Hue name
                boxpoints='all',  # Show all data points
                jitter=0.3,       # Spread points for visibility
                pointpos=0,       # Position points inside the box
                marker=dict(color=color_map[hue_val]),  # Color based on hue
                line=dict(width=line_width),  # Line styling
                fillcolor=color_map[hue_val],  # Fill color
                boxmean=True,  # Display the mean
                orientation='h'  # Horizontal box plots
            )
            
            traces_list.append(box_trace)
        
        return traces_list
    
    else:
        # Original functionality without hue
        unique_transcripts = long_expression_matrix[y].unique().to_list()
        
        # Sample names with transcript names
        long_expression_matrix = long_expression_matrix.with_columns(
            (pl.col(sample_id_column) + "_" + pl.col(y)).alias("trace_name")
        )
        
        # Initialize list of Box traces
        traces_list = []
    
        # Iterate over each unique transcript to create a Box trace
        for transcript in unique_transcripts:
            # Filter counts for the current transcript
            transcript_df = long_expression_matrix.filter(pl.col(y) == transcript)
            expression = transcript_df[x].to_list()
            sample_id = transcript_df[sample_id_column].to_list()
            
            # Create a Box trace
            box_trace = go.Box(
                y=expression,  # Numerical data for the box plot
                text=sample_id,  # Label for the trace
                name=transcript,
                pointpos=0,
                boxmean=True,  # Display the mean
                boxpoints='all',
                marker_color="black",  # Marker color
                line=dict(width=line_width),  # Line styling
                fillcolor="blue",
                orientation='h'  # Horizontal box plots
            )
            
            traces_list.append(box_trace)
        
        return traces_list
