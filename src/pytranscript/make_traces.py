import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
import polars as pl  # Imports Polars for data manipulation and analysis
from plotly.subplots import make_subplots  # Imports function for creating subplots in Plotly
from pytranscript.utils import check_df  # Imports utility function for data validation

def make_traces(
    data: pl.DataFrame,  # Input Polars dataframe containing transcript data
    x_start: str = "start",  # Column name for the start position (x-axis)
    x_end: str = "end",  # Column name for the end position (x-axis)
    y: str = "transcript_id",  # Column name for the transcript identifier (y-axis)
    strand: str = "strand",  # Column name for strand information (e.g., "+" or "-")
    type: str = "type",  # Column name indicating the type of feature (e.g., exon, intron)
    cds: str = "CDS",  # Value representing coding sequences (CDS)
    exon: str = "exon",  # Value representing exons
    intron: str = "intron",  # Value representing introns
    line_color: str = "black",  # Line color for the shapes (e.g., exon, CDS)
    fill_column: str = None,  # Optional column name for fill colors
    fill_color: str = "grey",  # Default fill color if no fill_column is specified
    intron_line_width: float = 0.5,  # Line width for intron shapes
    exon_line_width: float = 0.25,  # Line width for exon shapes
    opacity: float = 1,  # Opacity for all shapes
    arrow_color: str = "black",  # Arrow color for directional introns
    exon_height: float = 0.3,  # Height of exon shapes on the y-axis
    cds_height: float = 0.5,  # Height of CDS shapes on the y-axis
    arrow_height: float = 0.5,  # Height of directional arrows on the introns
    arrow_length: float = 1, ## Length of the arrow shapes (x-axis coordinate)
    arrow_line_width: float = 0.5  # Line width for the arrow shapes
) -> list:
   
    """
    Generates Plotly traces for visualizing transcript features such as exons, coding sequences (CDS), and introns.

    This function creates graphical shapes (rectangles for exons/CDS, lines for introns) to represent transcript structures 
    in a genomic visualization. It supports customizing the appearance of these features, including color, line width, 
    height, and the addition of directional arrows for introns.

    Parameters:
        data (pl.DataFrame): Polars DataFrame with transcript feature data.
        x_start (str, optional): Column name for feature start positions (default is 'start').
        x_end (str, optional): Column name for feature end positions (default is 'end').
        y (str, optional): Column name for transcript IDs to map to the y-axis (default is 'transcript_id').
        strand (str, optional): Column name for strand information (default is 'strand').
        type (str, optional): Column name indicating feature type (e.g., exon, intron, CDS) (default is 'type').
        cds (str, optional): Feature type representing coding sequences (default is 'CDS').
        exon (str, optional): Feature type representing exons (default is 'exon').
        intron (str, optional): Feature type representing introns (default is 'intron').
        line_color (str, optional): Line color for feature outlines (default is 'black').
        fill_column (str, optional): Column name for individual feature fill colors (optional).
        fill_color (str, optional): Default fill color for features if no fill_column is provided (default is 'grey').
        intron_line_width (float, optional): Line width for introns (default is 0.5).
        exon_line_width (float, optional): Line width for exons (default is 0.25).
        opacity (float, optional): Opacity of the shapes (default is 1.0).
        arrow_color (str, optional): Color of intron arrows (default is 'black').
        exon_height (float, optional): Vertical height of exon shapes (default is 0.3).
        cds_height (float, optional): Vertical height of CDS shapes (default is 0.5).
        arrow_height (float, optional): Vertical size of intron arrows (default is 0.5).
        arrow_length (float, optional): Length of intron arrows (default is 1.0).
        arrow_line_width (float, optional): Line width of intron arrows (default is 0.5).

    Returns:
        list: A list of Plotly trace objects for rendering genomic features in a plot.
    """

    # Validate required columns in the data
    if fill_column is None:
        check_df(data, [x_start, x_end, y, type, strand])
    else:
        check_df(data, [x_start, x_end, y, type, strand, fill_column])

    # Define trace lists for different types of features
    cds_traces = []  # Stores traces for CDS (coding sequences)
    intron_traces = []  # Stores traces for introns
    exon_traces = []  # Stores traces for exons

    # Get unique transcript IDs to map them to specific y-axis positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    
    # Create a dictionary that maps each unique transcript ID to a y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # Calculate the global maximum and minimum x-values (positions)
    global_max = max(
        data.select(pl.col(x_start).max()).item(),
        data.select(pl.col(x_end).max()).item()
    )

    global_min = min(
        data.select(pl.col(x_start).min()).item(),
        data.select(pl.col(x_end).min()).item()
    )

    # Calculate the total size of the x-axis range
    size = int(abs(global_max - global_min))

    # Handle the case where no specific fill column is provided
    if fill_column is None:
        fill = [fill_color] * len(data)  # Use the default fill color for all features
    else:
        fill = data[fill_column]  # Use the fill color from the specified column

    # Iterate over each row in the DataFrame to create traces for exons, CDS, and introns
    for idx, row in enumerate(data.iter_rows(named=True)):
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript

        # If the feature type is an exon, create a rectangle trace
        if row[type] == exon:
            trace = dict(
                type="rect",  # Rectangle trace for the exon
                x0=row[x_start],  # Start position on the x-axis
                x1=row[x_end],  # End position on the x-axis
                y0=y_pos - exon_height / 2,  # Bottom boundary of the rectangle on the y-axis
                y1=y_pos + exon_height / 2,  # Top boundary of the rectangle on the y-axis
                fillcolor=fill[idx],  # Fill color for the exon
                line=dict(color=line_color, width=exon_line_width),  # Border color and width
                opacity=opacity,  # Transparency level
            )
            exon_traces.append(trace)  # Append the trace to the exon list

        # If the feature type is a CDS, create a rectangle trace
        elif row[type] == cds:
            trace = dict(
                type="rect",  # Rectangle trace for the CDS
                x0=row[x_start],
                x1=row[x_end],
                y0=y_pos - cds_height / 2,  # Bottom boundary of the rectangle on the y-axis
                y1=y_pos + cds_height / 2,  # Top boundary of the rectangle on the y-axis
                fillcolor=fill[idx],  # Fill color for the CDS
                line=dict(color=line_color, width=exon_line_width),  # Border color and width
                opacity=opacity,
            )
            cds_traces.append(trace)  # Append the trace to the CDS list

        # If the feature type is an intron, create a line trace
        elif row[type] == intron:
            trace = dict(
                type="line",  # Line trace for the intron
                x0=row[x_start],
                x1=row[x_end],
                y0=y_pos,  # Introns are represented as horizontal lines
                y1=y_pos,
                line=dict(color=line_color, width=intron_line_width),  # Line color and width
                opacity=opacity,
            )
            intron_traces.append(trace)  # Append the trace to the intron list

            # Add directional arrows for introns if they are long enough
            if abs(row[x_start] - row[x_end]) > size / 25:
                arrow_x = (row[x_start] + row[x_end]) / 2  # Midpoint of the intron
                if row[strand] == "+":  # Positive strand (direction rightward)
                    arrow_trace = [
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x - (size / (150 / arrow_length)),  # Arrowhead to the left
                            y1=y_pos + arrow_height / 2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        ),
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x - (size / (150 / arrow_length)),  # Arrowhead to the left
                            y1=y_pos - arrow_height / 2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        )
                    ]
                elif row[strand] == "-":  # Negative strand (direction leftward)
                    arrow_trace = [
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x + (size / (150 / arrow_length)),  # Arrowhead to the right
                            y1=y_pos + arrow_height / 2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        ),
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x + (size / (150 / arrow_length)),  # Arrowhead to the right
                            y1=y_pos - arrow_height / 2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        )
                    ]

                # Append the arrow trace to the intron traces
                intron_traces.extend(arrow_trace)

    # Combine all traces (exons, CDS, introns, and arrows)
    traces = exon_traces + cds_traces + intron_traces

    return traces  # Return the list of traces

