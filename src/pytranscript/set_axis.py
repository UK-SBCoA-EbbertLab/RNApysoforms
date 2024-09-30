import polars as pl  # For data manipulation
import plotly.graph_objects as go  # For creating and updating Plotly figures
from typing import Union, List  # For type hinting
from pytranscript.utils import check_df  # Utility function for data validation

def set_axis(
    fig: go.Figure,  # Plotly figure object to be updated
    data: pl.DataFrame,  # Polars DataFrame containing the genomic data with 'start' and 'end' columns
    padding: int = 100,  # Optional padding to add space around the x-axis range, default is 100
    group_var: str = "transcript_id"  # Column name used to group data on the y-axis (default is 'transcript_id')
) -> go.Figure:
    """
    Adjust the x-axis range of the plot to align with genomic coordinates, and set the y-axis range
    based on the number of distinct groups (e.g., transcripts) in the data.
    
    Parameters:
    -----------
    fig : go.Figure
        The Plotly figure object where the data is being displayed.
    data : pl.DataFrame
        The DataFrame containing genomic data with 'start' and 'end' columns.
    padding : int, optional
        Extra space to add to the x-axis on both sides for better visualization. Default is 100.
    group_var : str, optional
        The column name used to group the data on the y-axis. Default is 'transcript_id'.
    
    Returns:
    --------
    go.Figure
        The updated figure object with corrected axes based on the genomic coordinates and grouping variable.
    
    Raises:
    -------
    ValueError
        If the 'start' and 'end' columns are not found in the data.
    """
    
    # Ensure the DataFrame contains the required columns ('start' and 'end')
    check_df(data, ["start", "end"])

    # Find the minimum start and maximum end values to define the x-axis range
    min_start = data['start'].min()  # Minimum value of the 'start' column
    max_end = data['end'].max()  # Maximum value of the 'end' column
    
    # Add padding to the genomic range for better visualization
    x_min = min_start - padding  # Minimum x-axis value with padding
    x_max = max_end + padding  # Maximum x-axis value with padding
    
    # Update the x-axis range in the Plotly figure based on the calculated genomic range
    fig.update_xaxes(range=[x_min, x_max])

    # Calculate the total number of distinct groups (e.g., transcripts) based on the group_var
    num_transcripts = data.n_unique(subset=group_var)

    # Update the y-axis range based on the number of distinct groups
    # Setting range to [-0.8, num_transcripts - 0.2] provides space for visual clarity
    fig.update_yaxes(range=[-0.8, (num_transcripts - 0.2)])

    return fig  # Return the updated figure
