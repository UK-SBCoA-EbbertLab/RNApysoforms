import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
import polars as pl  # Imports Polars for data manipulation and analysis

def geom_range(
    data,
    x_start="start",
    x_end="end",
    y="transcript_name",
    fill="grey",
    color="black",
    opacity=1,
    line_width=0.25,
    height=0.5,
):
    """
    Create a list of Plotly range (exon) trace objects.

    Parameters:
    -----------
    data : pl.DataFrame
        DataFrame containing the range (exon) data.
    x_start : str, optional
        Column name for the start position of exons. Default is 'start'.
    x_end : str, optional
        Column name for the end position of exons. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    fill : str or list, optional
        Fill color or a list of colors. Default is 'grey'.
    color : str, optional
        Color of the rectangle borders. Default is 'black'.
    opacity : float, optional
        Opacity of the rectangles. Default is 1.
    line_width : float, optional
        Width of the rectangle borders. Default is 0.25.
    height : float, optional
        Height of the rectangles. Default is 0.5.

    Returns:
    --------
    list
        A list of Plotly trace objects representing the exons.
    """
    traces = []  # Initialize an empty list to store Plotly trace objects (rectangles representing exons)

    # Get the unique values from the 'y' column (e.g., transcript names) to assign them distinct y-positions
    unique_y = data[y].unique(maintain_order=True).to_list()

    # Create a dictionary to map each unique transcript name to a numerical y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # If 'fill' is a single string, convert it into a list of the same color for all data points
    if isinstance(fill, str):
        fill = [fill] * len(data)

    # Iterate over each row in the DataFrame to create a Plotly rectangle for each exon
    for idx, row in enumerate(data.iter_rows(named=True)):
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript
        # Define the rectangle trace (an exon) with position and appearance attributes
        trace = dict(
            type="rect",  # Specifies that this trace is a rectangle
            x0=row[x_start],  # Start position of the exon (x-axis, left boundary)
            x1=row[x_end],  # End position of the exon (x-axis, right boundary)
            y0=y_pos - height / 2,  # Bottom boundary of the rectangle (y-axis)
            y1=y_pos + height / 2,  # Top boundary of the rectangle (y-axis)
            fillcolor=fill[idx],  # The color used to fill the rectangle (exon)
            line=dict(color=color, width=line_width),  # Border color and width of the rectangle
            opacity=opacity,  # Transparency of the rectangle
        )
        # Append the created trace to the 'traces' list
        traces.append(trace)

    # Return the list of Plotly trace objects
    return traces