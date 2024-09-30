import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
import polars as pl  # Imports Polars for data manipulation and analysis
from plotly.subplots import make_subplots


def make_traces(
    data,
    x_start="start",
    x_end="end",
    y="transcript_name",
    strand="strand",
    line_color="black",
    fill="grey",
    intron_line_width=0.5,
    exon_line_width=0.25,
    opacity=1,
    arrow_color="black",
    arrow_size=0.3,
    arrow_min_intron_length=100,
    exon_height=0.3,
    cds_height=0.5
):

    traces = []  # Initialize an empty list to store Plotly trace objects (rectangles representing exons/CDS)

    # Get the unique values from the 'y' column (e.g., transcript names) to assign them distinct y-positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    
    # Create a dictionary to map each unique transcript name to a numerical y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # If 'fill' is a single string, convert it into a list of the same color for all data points
    if isinstance(fill, str):
        fill = [fill] * len(exons)

    # Iterate over each row in the DataFrame to create a Plotly rectangle for each exon
    for idx, row in enumerate(data.iter_rows(named=True)):
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript

        if row["type"] == "exon":
            # Define the rectangle trace (an exon) with position and appearance attributes
            trace = dict(
                type="rect",  # Specifies that this trace is a rectangle
                x0=row[x_start],  # Start position of the exon (x-axis, left boundary)
                x1=row[x_end],  # End position of the exon (x-axis, right boundary)
                y0=y_pos - exon_height / 2,  # Bottom boundary of the rectangle (y-axis)
                y1=y_pos + exon_height / 2,  # Top boundary of the rectangle (y-axis)
                fillcolor=fill[idx],  # The color used to fill the rectangle (exon)
                line=dict(color=line_color, width=exon_line_width),  # Border color and width of the rectangle
                opacity=opacity,  # Transparency of the rectangle
            )

        elif row["type"] == "CDS":
            
            # Define the rectangle trace (an exon) with position and appearance attributes
            trace = dict(
                type="rect",  # Specifies that this trace is a rectangle
                x0=row[x_start],  # Start position of the exon (x-axis, left boundary)
                x1=row[x_end],  # End position of the exon (x-axis, right boundary)
                y0=y_pos - cds_height / 2,  # Bottom boundary of the rectangle (y-axis)
                y1=y_pos + cds_height / 2,  # Top boundary of the rectangle (y-axis)
                fillcolor=fill[idx],  # The color used to fill the rectangle (exon)
                line=dict(color=line_color, width=exon_line_width),  # Border color and width of the rectangle
                opacity=opacity,  # Transparency of the rectangle
            )

        elif row["type"] == "intron":
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
            
        # Append the created trace to the 'traces' list
        traces.append(trace)


