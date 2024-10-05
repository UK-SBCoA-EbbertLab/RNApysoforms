import plotly.graph_objects as go
import plotly.express as px
import polars as pl
from plotly.subplots import make_subplots
from RNApysoforms.utils import check_df
from typing import List


def make_traces(
    data: pl.DataFrame,
    x_start: str = "start",
    x_end: str = "end",
    y: str = "transcript_id",
    strand: str = "strand",
    feature_number: str= "exon_number",
    type: str = "type",
    cds: str = "CDS",
    exon: str = "exon",
    intron: str = "intron",
    line_color: str = "black",
    fill_color: str = "grey",
    hue: str = None,
    color_map: dict = None,
    color_palette: List[str] = px.colors.qualitative.Plotly,
    intron_line_width: float = 0.5,
    exon_line_width: float = 0.25,
    opacity: float = 1,
    arrow_color: str = "black",
    exon_height: float = 0.3,
    cds_height: float = 0.5,
    arrow_height: float = 0.5,
    arrow_length: float = 1,
    arrow_line_width: float = 0.5,
    is_hoverable: bool = True
) -> list:
    """
    Generates Plotly traces for visualizing transcript features like exons, introns, and coding sequences (CDS).

    This function creates graphical shapes (rectangles for exons and CDS, lines for introns) for transcript visualization 
    in genomic plots. It allows customization of appearance, such as color, line width, feature heights, and the option to 
    add directional arrows for introns to indicate strand direction. The function supports coloring features based on a 
    'hue' column, enabling grouping and coloring of features based on another variable.

    Parameters
    ----------
    data : pl.DataFrame
        A Polars DataFrame containing transcript feature data.
    x_start : str, optional
        Column name for feature start positions (x-axis), by default "start".
    x_end : str, optional
        Column name for feature end positions (x-axis), by default "end".
    y : str, optional
        Column name for transcript IDs to map to the y-axis, by default "transcript_id".
    strand : str, optional
        Column name for strand information (e.g., "+" or "-"), by default "strand".
    type : str, optional
        Column name indicating the feature type (e.g., exon, intron, CDS), by default "type".
    cds : str, optional
        Value representing coding sequences, by default "CDS".
    exon : str, optional
        Value representing exons, by default "exon".
    intron : str, optional
        Value representing introns, by default "intron".
    line_color : str, optional
        Line color for feature outlines, by default "black".
    fill_color : str, optional
        Default fill color if no `hue` is provided, by default "grey".
    hue : str, optional
        Column name for feature grouping to apply different fill colors.
    color_map : dict, optional
        A dictionary mapping unique hue values to specific colors. If None, a color map is generated using `color_palette`.
    color_palette : List[str], optional
        A list of colors to use when generating the color map if `color_map` is None and `hue` is provided, 
        by default `px.colors.qualitative.Plotly`.
    intron_line_width : float, optional
        Line width for introns, by default 0.5.
    exon_line_width : float, optional
        Line width for exons and CDS, by default 0.25.
    opacity : float, optional
        Opacity of all shapes, by default 1.
    arrow_color : str, optional
        Arrow color for directional introns, by default "black".
    exon_height : float, optional
        Vertical height of exon shapes, by default 0.3.
    cds_height : float, optional
        Vertical height of CDS shapes, by default 0.5.
    arrow_height : float, optional
        Vertical height of directional arrows on introns, by default 0.5.
    arrow_length : float, optional
        Length of intron arrows, by default 1.
    arrow_line_width : float, optional
        Line width for arrow shapes, by default 0.5.

    Returns
    -------
    list
        A list of Plotly trace objects representing exons, introns, CDS, and optional arrows.

    Examples
    --------
    Create traces for a genomic visualization:

    >>> import polars as pl
    >>> from RNApysoforms.utils import check_df
    >>> from RNApysoforms.plot import make_traces
    >>> df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx1"],
    ...     "start": [100, 400, 800],
    ...     "end": [200, 500, 900],
    ...     "type": ["exon", "intron", "CDS"],
    ...     "strand": ["+", "+", "+"],
    ...     "feature_group": ["group1", "group1", "group2"]
    ... })
    >>> traces = make_traces(df, hue="feature_group")
    >>> print(traces)
    
    This will return a list of Plotly traces that can be used to render the genomic features, 
    colored according to the 'feature_group' column.
    
    Notes
    -----
    - The function automatically adds arrows to introns if they are long enough and provides customization options 
      like arrow direction based on the strand column.
    - The y-axis positions of each feature are determined by the `transcript_id` column, and you can customize 
      colors and styles using the available parameters.
    - When a `hue` is provided, the features are colored according to the values in the `hue` column. If `color_map` 
      is not provided, it will generate a color map using `color_palette`.
    - The legend will display each unique hue value only once.
    - Arrows indicate strand direction: for positive strand ('+'), arrows point leftward; for negative strand ('-'), 
      arrows point rightward.
    """

    # Validate required columns in the data
    if hue is None:
        check_df(data, [x_start, x_end, y, type, strand, "exon_number"])
    else:
        check_df(data, [x_start, x_end, y, type, strand, hue])

    # Define trace lists for different types of features
    cds_traces = []    # Stores traces for CDS (coding sequences)
    intron_traces = []  # Stores traces for introns
    exon_traces = []   # Stores traces for exons

    # Get unique transcript IDs to map them to specific y-axis positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    
    # Create a dictionary that maps each unique transcript ID to a y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # Generate a color map if not provided and 'hue' is specified
    if ((color_map is None) and (hue is not None)):
        # Get unique values from the hue column to map to colors
        values_to_colormap = data[hue].unique(maintain_order=True).to_list()
        # Create a color map dictionary mapping hue values to colors from the palette
        color_map = {value: color for value, color in zip(values_to_colormap, color_palette)}
    # If 'hue' is not provided, set the color map to a single fill color
    elif hue is None:
        color_map = fill_color

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

    # Create a list to keep track of hue values already displayed in the legend
    displayed_hue_names = []

    # Iterate over each row in the DataFrame to create traces for exons, CDS, and introns
    for idx, row in enumerate(data.iter_rows(named=True)):
        
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript

        
        # Determine the fill color and legend name based on 'hue'
        if hue is None:
            exon_and_cds_color = fill_color
            hue_name = "Exon and/or CDS"
        else:
            exon_and_cds_color = color_map[row[hue]]
            hue_name = row[hue]

        # Determine whether to display the legend entry for this hue value
        if ((hue_name in displayed_hue_names) or (hue is None)):
            display_legend = False
        else: 
            display_legend = True
        
        if is_hoverable:
            ## Define hover template
            hovertemplate_text = (
                f"<b>Start:</b> {row[x_start]}<br>"
                f"<b>End:</b> {row[x_end]}<br>"
                f"<b>Feature Number:</b> {row[type]} {row[feature_number]}<br>"
                f"<b>{y}:</b> {row[y]}"
            )
        else:
            hovertemplate_text=None

        # If the feature type is an exon, create a scatter trace representing a rectangle
        if row[type] == exon:
            # Define the coordinates of the rectangle's corners
            x0 = row[x_start]  # Start position on the x-axis
            x1 = row[x_end]    # End position on the x-axis
            y0 = y_pos - exon_height / 2  # Bottom boundary of the rectangle on the y-axis
            y1 = y_pos + exon_height / 2  # Top boundary of the rectangle on the y-axis

            # Create the scatter trace
            trace = dict(
                type='scatter',
                mode='lines',
                x=[x0, x1, x1, x0, x0],  # Define the corners of the rectangle
                y=[y0, y0, y1, y1, y0],
                fill='toself',  # Fill the area enclosed by the lines
                fillcolor=exon_and_cds_color,  # Fill color for the exon
                line=dict(color=line_color, width=exon_line_width),  # Border color and width
                opacity=opacity,  # Transparency level
                name=hue_name,  # Legend name based on hue
                showlegend=display_legend,  # Show legend only once per hue
                hovertemplate=hovertemplate_text,  # Custom hover text
                hoverinfo='skip'  # Ensures the legend doesn't show up in the hover
            )
            exon_traces.append(trace)  # Append the trace to the exon list
            
            # Keep track of hue values that have been displayed in the legend
            if hue_name not in displayed_hue_names:
                displayed_hue_names.append(hue_name)

        # If the feature type is a CDS, create a scatter trace representing a rectangle
        elif row[type] == cds:
            # Define the coordinates of the rectangle's corners
            x0 = row[x_start]
            x1 = row[x_end]
            y0 = y_pos - cds_height / 2
            y1 = y_pos + cds_height / 2

            # Create the scatter trace
            trace = dict(
                type='scatter',
                mode='lines',
                x=[x0, x1, x1, x0, x0],  # Outline the rectangle
                y=[y0, y0, y1, y1, y0],
                fill='toself',
                fillcolor=exon_and_cds_color,
                line=dict(color=line_color, width=exon_line_width),
                opacity=opacity,
                name=hue_name,
                showlegend=display_legend,
                hovertemplate=hovertemplate_text,
                hoverinfo='skip'  # Ensures the legend doesn't show up in the hover
            )
            cds_traces.append(trace)  # Append the trace to the CDS list

            # Keep track of hue values that have been displayed in the legend
            if hue_name not in displayed_hue_names:
                displayed_hue_names.append(hue_name)

        # If the feature type is an intron, create a scatter trace representing a line
        elif row[type] == intron:
            trace = dict(
                type='scatter',
                mode='lines',
                x=[row[x_start], row[x_end]],  # Start and end positions on the x-axis
                y=[y_pos, y_pos],  # Constant y to create a horizontal line
                line=dict(color=line_color, width=intron_line_width),
                opacity=opacity,
                hovertemplate=hovertemplate_text,
                hoverinfo='skip',  # Ensures the legend doesn't show up in the hover
                showlegend=False  # Assuming introns don't need legend entries
            )
            intron_traces.append(trace)  # Append the trace to the intron list


            # Add directional arrows for introns if they are long enough
            if abs(row[x_start] - row[x_end]) > size / 25:
                arrow_x = (row[x_start] + row[x_end]) / 2  # Midpoint of the intron
                arrow_length_px = size / (150 / arrow_length)  # Calculate arrow length in plot units
                if row[strand] == "+":  # Positive strand (arrow pointing left)
                    x_coords = [
                        arrow_x, arrow_x - arrow_length_px, None,
                        arrow_x, arrow_x - arrow_length_px
                    ]
                    y_coords = [
                        y_pos, y_pos + arrow_height / 2, None,
                        y_pos, y_pos - arrow_height / 2
                    ]
                elif row[strand] == "-":  # Negative strand (arrow pointing right)
                    x_coords = [
                        arrow_x, arrow_x + arrow_length_px, None,
                        arrow_x, arrow_x + arrow_length_px
                    ]
                    y_coords = [
                        y_pos, y_pos + arrow_height / 2, None,
                        y_pos, y_pos - arrow_height / 2
                    ]

                arrow_trace = go.Scatter(
                    x=x_coords,
                    y=y_coords,
                    mode='lines',
                    line=dict(color=arrow_color, width=arrow_line_width),
                    opacity=opacity,
                    showlegend=False,
                    hoverinfo="skip"
                )

                # Append the arrow trace to the intron traces
                intron_traces.append(arrow_trace)

    # Combine all traces (exons, CDS, introns, and arrows)
    traces = exon_traces + cds_traces + intron_traces

    return traces  # Return the list of traces
