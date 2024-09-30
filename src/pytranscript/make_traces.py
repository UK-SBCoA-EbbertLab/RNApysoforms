import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
import polars as pl  # Imports Polars for data manipulation and analysis
from plotly.subplots import make_subplots


def make_traces(
    data, ## Polars dataframe
    x_start="start",
    x_end="end",
    y="transcript_name",
    type="type",
    cds="CDS", 
    exon="exon",
    intron="intron",
    strand="strand",
    line_color="black",
    fill_column=None,
    fill_color="grey",
    intron_line_width=0.5,
    exon_line_width=0.25,
    opacity=1,
    arrow_color="black",
    exon_height=0.3,
    cds_height=0.5,
    arrow_height=0.5,
    arrow_line_width=0.5
):
    
    ## Order data starting with intron, then CDS, then exon so traces are made in correct order
    data = data.with_columns(
        pl.when(pl.col(type) == exon).then(0)
        .when(pl.col(type) == cds).then(1)
        .when(pl.col(type) == intron).then(2)
        .otherwise(3).alias("type_order")
    ).sort("type_order").drop("type_order")

    ## Define trace lists
    cds_traces = []
    intron_traces = []
    exon_traces = []

    # Get the unique values from the 'y' column (e.g., transcript names) to assign them distinct y-positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    
    # Create a dictionary to map each unique transcript name to a numerical y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    
    # Calculate global maximum
    global_max = max(
        data.select(pl.col(x_start).max()).item(),
        data.select(pl.col(x_end).max()).item()
    )

    # Calculate global minimum
    global_min = min(
        data.select(pl.col(x_start).min()).item(),
        data.select(pl.col(x_end).min()).item()
    )

    # Calculate size
    size = int(abs(global_max - global_min))

    # If 'fill' is a single string, convert it into a list of the same color for all data points
    if fill_column == None:
        fill=[fill_color] * len(data)
    else:
        fill = data[fill_column]

    # Iterate over each row in the DataFrame to create a Plotly rectangle for each exon
    for idx, row in enumerate(data.iter_rows(named=True)):
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript

        if row[type] == exon:
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
            exon_traces.append(trace)

        elif row[type] == cds:
            
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
            cds_traces.append(trace)

        elif row[type] == intron:
            # Create a line trace for the intron
            trace = dict(
                type="line",
                x0=row[x_start],
                x1=row[x_end],
                y0=y_pos,
                y1=y_pos,
                line=dict(color=line_color, width=intron_line_width),
                opacity=opacity,
            )
            intron_traces.append(trace)


            # Ensure arrows are not too close to the intron ends
            if (
                abs(row[x_start] - row[x_end]) > size/25
            ):
                 
                arrow_x = (row[x_start] + row[x_end])/2
                if row[strand] == "+":
                    # Create the two line shapes for the arrow (like `>` shape)
                    arrow_trace = [
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x - size/250,
                            y1=y_pos + arrow_height/2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        ),
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x - size/250,
                            y1=y_pos - arrow_height/2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        )
                    ]

                elif row[strand] == "-":
                    # Create the two line shapes for the arrow (like `>` shape)
                    arrow_trace = [
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x + size/250,
                            y1=y_pos + arrow_height/2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        ),
                        dict(
                            type="line",
                            x0=arrow_x,
                            y0=y_pos,
                            x1=arrow_x + size/250,
                            y1=y_pos - arrow_height/2,
                            line=dict(color=arrow_color, width=arrow_line_width),
                            opacity=opacity,
                        )
                    ]


                intron_traces.extend(arrow_trace)
        
        traces = exon_traces + cds_traces + intron_traces
            
    return traces


