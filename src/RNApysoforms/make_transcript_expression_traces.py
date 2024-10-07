import plotly.graph_objects as go
import polars as pl
from typing import List

def make_transcript_expression_traces(
    long_expression_matrix: pl.DataFrame,
    y: str = "transcript_id",
    x: str = "counts",
    line_width: float = 0.5
) -> List[go.Box]:
    """
    Create a list of Plotly Box traces for each unique transcript.

    Parameters:
    - long_expression_matrix (pl.DataFrame): DataFrame containing transcript expression data.
    - y (str): Column name for transcript IDs (categorical axis).
    - x (str): Column name for counts (numerical axis).

    Returns:
    - List[go.Box]: List of Plotly Box traces.
    """
    
    # Ensure the necessary columns exist
    if y not in long_expression_matrix.columns or x not in long_expression_matrix.columns:
        raise ValueError(f"Columns '{y}' and/or '{x}' not found in the DataFrame.")
    
    # Create a consistent y-axis mapping
    unique_transcripts = long_expression_matrix[y].unique(maintain_order=True).to_list()
    y_dict = {val: i for i, val in enumerate(unique_transcripts)}
    
    
    # Initialize list of Box traces
    traces_list = []

    # Iterate over each unique transcript to create a Box trace
    for transcript in unique_transcripts:

        y_pos = y_dict[transcript]  # Get the corresponding y-position for the current transcript

        # Filter counts for the current transcript
        expression = long_expression_matrix.filter(pl.col(y) == transcript)[x].to_list()

        ## Sample Names
        
        # Create a Box trace
        box_trace = go.Box(
            x=expression,  # Numerical data for the box plot
            #y=[y_pos] * len(expression), 
            name=transcript,  # Label for the trace
            boxmean=True,  # Display the mean
            boxpoints='all',
            marker_color="black",  # Marker color
            line=dict(width=line_width),  # Line styling
            fillcolor="blue"
        )

        print(box_trace)
        
        traces_list.append(box_trace)
    
    return traces_list
