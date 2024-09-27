import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def geom_intron(
    data,
    x_start="start",
    x_end="end",
    y="transcript_name",
    strand="strand",
    color="black",
    line_width=0.5,
    opacity=1,
    arrow_color="black",
    arrow_size=0.3,
    arrow_min_intron_length=100,
):
    """
    Create a list of Plotly intron trace objects with strand arrows.

    Parameters:
    -----------
    data : pl.DataFrame
        DataFrame containing the intron data.
    x_start : str, optional
        Column name for the start position of introns. Default is 'start'.
    x_end : str, optional
        Column name for the end position of introns. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    strand : str, optional
        Column name for the strand information. Default is 'strand'.
    color : str, optional
        Color of the intron lines. Default is 'black'.
    line_width : float, optional
        Width of the intron lines. Default is 0.5.
    opacity : float, optional
        Opacity of the intron lines. Default is 1.
    arrow_color : str, optional
        Color of the strand arrows. Default is 'black'.
    arrow_size : float, optional
        Size of the strand arrows. Default is 0.3.
    arrow_min_intron_length : int, optional
        Minimum intron length to draw an arrow. Default is 100.

    Returns:
    --------
    list
        A list of Plotly trace objects representing the introns with optional strand arrows.
    """
    traces = []

    # Get unique values of the y-axis column and create a mapping to positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # Iterate over each row in the DataFrame
    for row in data.iter_rows(named=True):
        y_pos = y_dict[row[y]]
        
        # Create a line trace for the intron
        trace = dict(
            type="line",
            x0=row[x_start],
            x1=row[x_end],
            y0=y_pos,
            y1=y_pos,
            line=dict(color=color, width=line_width),
            opacity=opacity,
        )
        traces.append(trace)

        # Define arrow direction based on strand
        arrow_direction = ">" if row[strand] == "+" else "<"

        # Check if intron length exceeds minimum length for drawing arrows
        intron_length = row[x_end] - row[x_start]
        if intron_length > arrow_min_intron_length:
            if arrow_min_intron_length < 1:
                arrow_min_intron_length = 100

            num_arrows = int(intron_length / arrow_min_intron_length)
            for i in range(num_arrows):
                arrow_x = row[x_start] + i * (intron_length / num_arrows)
                # Ensure arrows are not too close to the intron ends
                if (
                    abs(arrow_x - row[x_end]) > 100
                    and abs(arrow_x - row[x_start]) > 100
                ):
                    # Create a scatter trace for the arrow
                    arrow_trace = go.Scatter(
                        x=[arrow_x],
                        y=[y_pos],
                        mode="markers",
                        marker=dict(
                            symbol="arrow-right"
                            if row[strand] == "+"
                            else "arrow-left",
                            size=(arrow_size * 10),
                            color=arrow_color,
                        ),
                        showlegend=False,
                        hoverinfo="none",
                    )
                    traces.append(arrow_trace)

    return traces