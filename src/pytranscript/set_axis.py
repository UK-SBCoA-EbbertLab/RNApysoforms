import polars as pl
import plotly.graph_objects as go
from typing import Union, List
from rna-pysoforms.utils import check_df

def set_axis(
    fig: go.Figure,
    data: pl.DataFrame,
    padding: int = 100,
    group_var: str = "transcript_id"
) -> go.Figure:
    
    """
    Adjusts the x-axis and y-axis ranges of a Plotly figure to fit genomic coordinates and groupings.

    This function updates the x-axis range based on the genomic coordinates from the `start` and `end` columns 
    of the input DataFrame, applying optional padding for better visualization. It also configures the y-axis 
    based on the number of unique transcript groups (or other specified groups) to ensure clear, readable plots.

    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to be updated.
    data : pl.DataFrame
        A Polars DataFrame containing genomic data, which must include 'start' and 'end' columns.
    padding : int, optional
        Amount of padding to add around the x-axis range for visualization purposes (default is 100).
    group_var : str, optional
        The column name used to group data for the y-axis (default is 'transcript_id').

    Returns
    -------
    go.Figure
        The updated Plotly figure with adjusted x-axis and y-axis ranges.

    Raises
    ------
    ValueError
        If the 'start' or 'end' columns are missing from the input DataFrame.

    Examples
    --------
    Adjust the axes of a Plotly figure based on genomic data:

    >>> import polars as pl
    >>> import plotly.graph_objects as go
    >>> from rna-pysoforms.plot import set_axis
    >>> df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
    ...     "start": [100, 200, 300, 400],
    ...     "end": [150, 250, 350, 450]
    ... })
    >>> fig = go.Figure()
    >>> updated_fig = set_axis(fig, df)
    >>> updated_fig.show()
    
    This will adjust the x-axis to fit the genomic ranges and set the y-axis according to the unique transcript IDs.
    
    Notes
    -----
    - The x-axis range is defined by the minimum `start` and maximum `end` values, with optional padding applied.
    - The y-axis range is calculated based on the number of unique groups defined by `group_var`, such as transcripts.
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
