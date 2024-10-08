import plotly.graph_objects as go
import polars as pl
from typing import List, Optional, Dict
from RNApysoforms.utils import check_df
import plotly.express as px

def make_transcript_expression_traces(
    expression_matrix: pl.DataFrame,
    y: str = "transcript_id",
    expression_columns: List[str] = ["counts"],
    sample_id_column: str = "sample_id",
    hue: Optional[str] = None,  # New hue parameter
    fill_color="grey",
    line_width: float = 0.5,
    color_palette: List[str] = px.colors.qualitative.Plotly_r,
    color_map: dict = None,  # Optional color map for hues    
) -> List[go.Box]:
    """
    Create a list of Plotly Box traces for each unique transcript, colored and dodged by a hue variable.

    Parameters:
    - expression_matrix (pl.DataFrame): DataFrame containing transcript expression data.
    - y (str): Column name for transcript IDs (categorical axis).
    - x (str): Column name for counts (numerical axis).
    - hue (str, optional): Column name for the hue (grouping) variable.
    - line_width (float): Width of the box plot lines.
    - sample_id_column (str): Column name for sample IDs.
    - color_map (dict, optional): Dictionary mapping hue values to colors.

    Returns:
    - List[go.Box]: List of Plotly Box traces.
    """

    ## TODO: Validate expression_columns as a list

    # Ensure the necessary columns exist
    required_columns = [y, sample_id_column] + expression_columns
    if hue is not None:
        required_columns.append(hue)
    check_df(expression_matrix, required_columns)
    
    # Generate a color map if not provided and 'hue' is specified
    if (color_map is None) and (hue is not None):
        # Get unique values from the hue column to map to colors
        values_to_colormap = expression_matrix[hue].unique(maintain_order=True).to_list()
        # Create a color map dictionary mapping hue values to colors from the palette
        color_map = {value: color for value, color in zip(values_to_colormap, color_palette)}
    elif hue is None:
        # If 'hue' is not provided, use a single fill color for all features
        color_map = fill_color

    ## 
    show_legend = True

    if hue is not None:

        ## List all unique hue values
        unique_hues = expression_matrix.get_column(hue).unique()
                
        # Initialize list of Box traces
        traces_list = []
        
        for x in expression_columns:

            x_traces_list = []

            # Iterate over each unique hue to create Box traces
            for hue_val in unique_hues:

                # Filter data for the current hue
                hue_filtered_df = expression_matrix.filter(pl.col(hue) == hue_val)
                
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
                    orientation='h',  # Horizontal box plots
                    legendgroup=str(hue_val),
                    showlegend=show_legend
                )
                
                x_traces_list.append(box_trace)
            
            show_legend=False    
            traces_list.append(x_traces_list)
            
        return traces_list
        
    else:
        
        # Original functionality without hue
        unique_transcripts = expression_matrix[y].unique().to_list()
        
        # Sample names with transcript names
        expression_matrix = expression_matrix.with_columns(
            (pl.col(sample_id_column) + "_" + pl.col(y)).alias("trace_name")
        )
        
        # Initialize list of Box traces
        traces_list = []

        for x in expression_columns:

            x_traces_list = []

            # Iterate over each unique transcript to create a Box trace
            for transcript in unique_transcripts:
                # Filter counts for the current transcript
                transcript_df = expression_matrix.filter(pl.col(y) == transcript)
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
                    fillcolor=fill_color,
                    orientation='h',  # Horizontal box plots
                    showlegend=show_legend
                )
                
                x_traces_list.append(box_trace)
            
            show_legend = False
            traces_list.append(x_traces_list)
        
        return traces_list
