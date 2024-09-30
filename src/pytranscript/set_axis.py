import polars as pl
import plotly.graph_objects as go
from typing import Union, List
from pytranscript.utils import check_df ## Import utils

def set_axis(fig, data: pl.DataFrame, padding=100, group_var="transcript_id"):
    """
    Adjust the x-axis range of the plot to align with genomic coordinates.

    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The plotly figure object where the data is being displayed.
    data : pl.DataFrame
        The DataFrame containing exon information with 'start' and 'end' columns.
    padding : int, optional
        Extra space to add to the x-axis on both sides for better visualization. Default is 100.

    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The updated figure object with corrected axes.
    """
    ## Check required columns:
    check_df(data, ["start", "end"])

    # Find the minimum start and maximum end values from both exon and intron data
    min_start = data['start'].min()
    max_end = data['end'].max()
    
    # Add padding to provide space around the genomic coordinates
    x_min = min_start - padding
    x_max = max_end + padding
    
    # Update the x-axis range of the plot to reflect the correct genomic range
    fig.update_xaxes(range=[x_min, x_max])

    # Total number of transcripts
    num_transcripts = data.n_unique(subset=group_var)

    # Update the y-axis range of the plot to reflect the number of transcripts
    fig.update_yaxes(range=[-0.8, (num_transcripts - 0.2)])

    return fig
