import polars as pl
import plotly.graph_objects as go
from typing import Union, List
from RNApysoforms.utils import check_df

def set_axis(
    fig: go.Figure,
    data: pl.DataFrame,
    x_start: str = "start",
    x_end: str = "end",
    padding: int = 100,
    transcript_id_column: str = "transcript_id"
) -> go.Figure:
    """
    Adjusts the x-axis and y-axis ranges of a Plotly figure to fit genomic coordinates and groupings.

    This function updates the x-axis range based on the genomic coordinates from the columns specified by `x_start` and `x_end` in the input DataFrame,
    applying optional padding for better visualization. It also configures the y-axis based on the number of unique groups defined by `transcript_id_column`
    (e.g., transcript IDs) to ensure clear and readable plots.

    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to be updated.
    data : pl.DataFrame
        A Polars DataFrame containing genomic data, which must include the columns specified by `x_start` and `x_end`.
    x_start : str, optional
        Column name representing the start position of genomic features, by default "start".
    x_end : str, optional
        Column name representing the end position of genomic features, by default "end".
    padding : int, optional
        Amount of padding to add around the x-axis range for visualization purposes (default is 100).
    transcript_id_column : str, optional
        The column name used to group data for the y-axis (e.g., transcript IDs), by default "transcript_id".

    Returns
    -------
    go.Figure
        The updated Plotly figure with adjusted x-axis and y-axis ranges.

    Raises
    ------
    TypeError
        - If `fig` is not an instance of `go.Figure`.
        - If `data` is not a Polars DataFrame.
    ValueError
        - If the columns specified by `x_start` or `x_end` are missing from the input DataFrame.

    Examples
    --------
    Adjust the axes of a Plotly figure based on genomic data:

    >>> import polars as pl
    >>> import plotly.graph_objects as go
    >>> from RNApysoforms.plot import set_axis
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
    - The x-axis range is defined by the minimum value in the `x_start` column and the maximum value in the `x_end` column, with optional padding applied.
    - The y-axis range is calculated based on the number of unique groups defined by `transcript_id_column`, such as transcripts, allowing multiple groups to be visualized clearly.
    - Customizing `x_start` and `x_end` allows flexibility in handling different genomic coordinate columns.
    - Ensure that the columns specified by `x_start` and `x_end` contain numeric values representing genomic positions.
    - Padding enhances visualization by preventing features from being too close to the plot boundaries.
    - The function does not modify the figure's layout beyond updating the axis ranges.

    """
    
    # Check if fig is an instance of go.Figure
    if not isinstance(fig, go.Figure):
        raise TypeError(f"Expected fig to be of type go.Figure, got {type(fig)}")
    
    # Check if data is a Polars DataFrame
    if not isinstance(data, pl.DataFrame):
        raise TypeError(
            f"Expected data to be of type pl.DataFrame, got {type(data)}" +
            "\n You can use polars_df = pl.from_pandas(pandas_df) to convert a pandas df into a polars df"
        )
    
    # Ensure the DataFrame contains the required columns specified by x_start and x_end
    check_df(data, [x_start, x_end])
    
    # Find the minimum start and maximum end values to define the x-axis range
    min_start = data[x_start].min()  # Minimum value of the `x_start` column
    max_end = data[x_end].max()      # Maximum value of the `x_end` column
    
    # Add padding to the genomic range for better visualization
    x_min = min_start - padding  # Minimum x-axis value with padding
    x_max = max_end + padding    # Maximum x-axis value with padding
    
    # Update the x-axis range in the Plotly figure based on the calculated genomic range
    fig.update_xaxes(range=[x_min, x_max])
    
    # Calculate the total number of distinct groups (e.g., transcripts) based on the transcript_id_column
    num_transcripts = data[transcript_id_column].n_unique()
    
    # Update the y-axis range based on the number of distinct groups
    # Setting range to [-0.8, num_transcripts - 0.2] provides space for visual clarity
    fig.update_yaxes(range=[-0.8, (num_transcripts - 0.2)])
    
    return fig  # Return the updated figure
