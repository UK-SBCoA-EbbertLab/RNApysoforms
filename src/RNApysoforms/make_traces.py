import plotly.graph_objects as go
import plotly.express as px
import polars as pl
from plotly.subplots import make_subplots
from RNApysoforms.shorten_gaps import shorten_gaps
from RNApysoforms.utils import check_df
from typing import List

def make_traces(
    data: pl.DataFrame,
    y: str = "transcript_id",
    x_start: str = "start",
    x_end: str = "end",
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
    exon_height: float = 0.3,
    cds_height: float = 0.5,
    arrow_height: float = 0.5,
    arrow_length: float = 1,
    is_hoverable: bool = True,
    hover_start: str = "start",
    hover_end: str = "end"
) -> list:
    """
    Generates Plotly traces for visualizing transcript features such as exons, introns, and coding sequences (CDS).

    This function creates graphical shapes—rectangles for exons and CDS, and lines for introns—for transcript visualization 
    in genomic plots. It allows extensive customization of appearance, including colors, line widths, feature heights, 
    and the option to add directional arrows to introns to indicate strand direction. Features can be colored based on a 
    specified 'hue' column, facilitating grouping and differentiation based on another variable. Hover information can 
    display details like feature type, feature number, and start and end positions. Additionally, the function supports 
    rescaling genomic coordinates to shorten long gaps, enhancing visualization clarity.

    **Required Columns in `data`:**
    - Columns specified by `y`, `x_start`, `x_end`.
    - `type`: Indicates the feature type (e.g., exon, intron, CDS).
    - `strand`: Indicates the strand direction ('+' or '-').
    - `exon_number`: Numerical identifier for exons.
    - `seqnames`: Chromosome or sequence name.
    - If `hue` is provided, the column specified by `hue` must also be present.

    Parameters
    ----------
    data : pl.DataFrame
        A Polars DataFrame containing transcript feature data. Must include columns for transcript ID, start and end 
        positions, feature type, strand direction, exon number, and sequence names. If `hue` is specified, the 
        corresponding column must also be present.
    y : str, optional
        Column name for transcript IDs to map to the y-axis, by default "transcript_id".
    x_start : str, optional
        Column name representing the start position of features, by default "start".
    x_end : str, optional
        Column name representing the end position of features, by default "end".
    cds : str, optional
        Value in the 'type' column representing coding sequences, by default "CDS".
    exon : str, optional
        Value in the 'type' column representing exons, by default "exon".
    intron : str, optional
        Value in the 'type' column representing introns, by default "intron".
    line_color : str, optional
        Color for the outlines of all features, by default "black".
    fill_color : str, optional
        Default fill color for features if no `hue` is provided, by default "grey".
    hue : str, optional
        Column name for feature grouping to apply different fill colors. If provided, features will be colored 
        based on the unique values in this column, enabling visual differentiation of groups.
    color_map : dict, optional
        A dictionary mapping unique `hue` values to specific colors. If `None`, a color map is generated using 
        `color_palette`.
    color_palette : List[str], optional
        A list of colors to use when generating the color map if `color_map` is `None` and `hue` is provided. 
        Defaults to `px.colors.qualitative.Plotly`.
    intron_line_width : float, optional
        Line width for intron traces, by default 0.5.
    exon_line_width : float, optional
        Line width for exon and CDS traces, by default 0.25.
    opacity : float, optional
        Opacity level for all feature shapes, ranging from 0 (transparent) to 1 (opaque), by default 1.
    exon_height : float, optional
        Vertical height of exon rectangles, by default 0.3.
    cds_height : float, optional
        Vertical height of CDS rectangles, by default 0.5.
    arrow_height : float, optional
        Vertical height of directional arrows on introns, by default 0.5.
    arrow_length : float, optional
        Length of the directional arrows on introns, by default 1.
    is_hoverable : bool, optional
        If `True`, features will display hover information including feature type, feature number, start and 
        end positions; by default `True`.
    hover_start : str, optional
        Column name representing the start position for hover information, by default "start".
    hover_end : str, optional
        Column name representing the end position for hover information, by default "end".

    Returns
    -------
    list
        A list of Plotly trace objects representing exons, introns, CDS, and optional arrows for visualization.

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
    ...     "feature_group": ["group1", "group1", "group2"],
    ...     "exon_number": [1, 1, 2],
    ...     "seqnames": ["chr1", "chr1", "chr1"]
    ... })
    >>> traces = make_traces(df, hue="feature_group")
    >>> print(traces)
        
    This will return a list of Plotly traces that can be used to render the genomic features, 
    colored according to the 'feature_group' column.

    Notes
    -----
    - The function automatically adds arrows to introns if they are sufficiently long, with arrow direction 
      based on the `strand` column.
    - Y-axis positions of each feature are determined by the `y` column, allowing multiple transcripts to be 
      visualized on the same plot.
    - When a `hue` is provided, features are colored according to the values in the `hue` column. If `color_map` 
      is not provided, a color map is generated using `color_palette`.
    - The legend will display each unique `hue` value only once to avoid redundancy.
    - Arrows indicate strand direction: for positive strand ('+'), arrows point leftward; for negative strand ('-'), 
      arrows point rightward.
    - The `exon_number` column is required when `hue` is `None` and is used to display feature numbers in hover 
      information.
    - The `seqnames` column is required for displaying chromosome or sequence names in hover information.
    - The `is_hoverable` parameter controls whether hover information is displayed for features. When `is_hoverable` 
      is `True`, features will display hover information including feature type, feature number, start and end positions.

    Raises
    ------
    TypeError
        If the `data` parameter is not a Polars DataFrame.
    ValueError
        If required columns are missing from the `data` DataFrame based on the provided parameters.
    """
    
    # Check if data is a Polars DataFrame
    if not isinstance(data, pl.DataFrame):
        raise TypeError(
            f"Expected data to be of type pl.DataFrame, got {type(data)}. "
            "You can use polars_df = pl.from_pandas(pandas_df) to convert a pandas DataFrame into a Polars DataFrame."
        )
    
    # Validate required columns in the data based on the presence of 'hue'
    if hue is None:
        # Required columns when 'hue' is not specified
        required_columns = [x_start, x_end, y, "type", "strand", "exon_number", "seqnames"]
    else:
        # Required columns when 'hue' is specified
        required_columns = [x_start, x_end, y, "type", "strand", "exon_number", hue, "seqnames"]
    
    check_df(data, required_columns)
    
    # Initialize lists to store traces for different feature types
    cds_traces = []      # Stores traces for CDS (coding sequences)
    intron_traces = []   # Stores traces for introns
    exon_traces = []     # Stores traces for exons
    
    # Get unique transcript IDs to map them to specific y-axis positions
    unique_y = data[y].unique(maintain_order=True).to_list()
    
    # Create a dictionary that maps each unique transcript ID to a y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}
    
    # Generate a color map if not provided and 'hue' is specified
    if (color_map is None) and (hue is not None):
        # Get unique values from the hue column to map to colors
        values_to_colormap = data[hue].unique(maintain_order=True).to_list()
        # Create a color map dictionary mapping hue values to colors from the palette
        color_map = {value: color for value, color in zip(values_to_colormap, color_palette)}
    elif hue is None:
        # If 'hue' is not provided, use a single fill color for all features
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
    for row in data.iter_rows(named=True):
        
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript
    
        # Determine the fill color and legend name based on 'hue'
        if hue is None:
            exon_and_cds_color = fill_color
            hue_name = "Exon and/or CDS"
        else:
            exon_and_cds_color = color_map.get(row[hue], fill_color)
            hue_name = row[hue]
    
        # Determine whether to display the legend entry for this hue value
        if (hue_name in displayed_hue_names) or (hue is None):
            display_legend = False
        else: 
            display_legend = True
        
        if is_hoverable:
            # Define hover template with feature type, number, start, and end positions for each row
            feature_size = abs(row[hover_end] - row[hover_start])
            hovertemplate_text = (
                f"<b>Feature Type:</b> {row['type']}<br>"
                f"<b>Feature Number:</b> {row.get('exon_number', 'N/A')}<br>"
                f"<b>Chromosome:</b> {row['seqnames']}<br>"
                f"<b>Start:</b> {row[hover_start]}<br>"
                f"<b>End:</b> {row[hover_end]}<br>"
                f"<b>Size:</b> {feature_size}<br>"
                "<extra></extra>"
            )
        else:
            hovertemplate_text = None
    
        # Create trace based on the feature type
        if row["type"] == exon:
            # Define the coordinates of the rectangle's corners for exon
            x0 = row[x_start]                   # Start position on the x-axis
            x1 = row[x_end]                     # End position on the x-axis
            y0 = y_pos - exon_height / 2        # Bottom boundary of the rectangle on the y-axis
            y1 = y_pos + exon_height / 2        # Top boundary of the rectangle on the y-axis
    
            # Create the scatter trace for the exon
            trace = dict(
                type='scatter',
                mode='lines+markers',
                x=[x0, x1, x1, x0, x0],          # Define the corners of the rectangle
                y=[y0, y0, y1, y1, y0],
                fill='toself',                    # Fill the area enclosed by the lines
                fillcolor=exon_and_cds_color,     # Fill color for the exon
                line=dict(color=line_color, width=exon_line_width),  # Border color and width
                marker=dict(opacity=0),           # Hide markers
                opacity=opacity,                  # Transparency level
                name=hue_name,                    # Legend name based on hue
                legendgroup=hue_name,
                showlegend=display_legend,        # Show legend only once per hue
                hovertemplate=hovertemplate_text, # Custom hover text
                hoverlabel=dict(namelength=-1),
                hoveron='fills+points'
            )
            exon_traces.append(trace)  # Append the trace to the exon list
            
            # Keep track of hue values that have been displayed in the legend
            if hue_name not in displayed_hue_names:
                displayed_hue_names.append(hue_name)
    
        elif row["type"] == cds:
            # Define the coordinates of the rectangle's corners for CDS
            x0 = row[x_start]
            x1 = row[x_end]
            y0 = y_pos - cds_height / 2
            y1 = y_pos + cds_height / 2
    
            # Create the scatter trace for the CDS
            trace = dict(
                type='scatter',
                mode='lines+markers',
                x=[x0, x1, x1, x0, x0],          # Outline the rectangle
                y=[y0, y0, y1, y1, y0],
                fill='toself',
                fillcolor=exon_and_cds_color,
                line=dict(color=line_color, width=exon_line_width),
                marker=dict(opacity=0),           # Hide markers
                opacity=opacity,
                name=hue_name,
                legendgroup=hue_name,
                showlegend=display_legend,
                hovertemplate=hovertemplate_text, # Custom hover text
                hoverlabel=dict(namelength=-1),
                hoveron='fills+points'
            )
            cds_traces.append(trace)  # Append the trace to the CDS list
    
            # Keep track of hue values that have been displayed in the legend
            if hue_name not in displayed_hue_names:
                displayed_hue_names.append(hue_name)
                
        elif row["type"] == intron:
            # Define the coordinates for the intron line
            x_intron = [row[x_start], row[x_end]]  # Start and end positions on the x-axis
            y_intron = [y_pos, y_pos]              # Constant y to create a horizontal line
    
            # Add directional arrows for introns if they are long enough
            if abs(row[x_start] - row[x_end]) > size / 25:
                arrow_x = (row[x_start] + row[x_end]) / 2  # Midpoint of the intron
                # Calculate arrow length in data units
                arrow_length_px = size / (150 / arrow_length) if arrow_length != 0 else 0  
    
                if row["strand"] == "+":  # Positive strand (arrow pointing left)
                    x_arrow = [
                        None,                              # Break before starting the arrow
                        arrow_x, arrow_x - arrow_length_px, None,
                        arrow_x, arrow_x - arrow_length_px, None  # Break after arrow
                    ]
                    y_arrow = [
                        None,
                        y_pos, y_pos + arrow_height / 2, None,
                        y_pos, y_pos - arrow_height / 2, None
                    ]
                
                elif row["strand"] == "-":  # Negative strand (arrow pointing right)
                    x_arrow = [
                        None,
                        arrow_x, arrow_x + arrow_length_px, None,
                        arrow_x, arrow_x + arrow_length_px, None
                    ]
                    y_arrow = [
                        None,
                        y_pos, y_pos + arrow_height / 2, None,
                        y_pos, y_pos - arrow_height / 2, None
                    ]
    
                # Combine the full intron line with None and arrow coordinates
                x_combined = x_intron + x_arrow
                y_combined = y_intron + y_arrow
            else:
                # If no arrow is added, use just the full intron coordinates
                x_combined = x_intron
                y_combined = y_intron
    
            # Create the scatter trace for the intron with optional arrows
            trace = dict(
                type='scatter',
                mode='lines',
                x=x_combined,
                y=y_combined,
                line=dict(color=line_color, width=intron_line_width),
                opacity=opacity,
                hovertemplate=hovertemplate_text,  # Make intron hoverable with this
                showlegend=False                   # Introns don't need legend entries
            )
            intron_traces.append(trace)  # Append the trace to the intron list
    
    # Combine all traces (exons, CDS, introns)
    traces = exon_traces + cds_traces + intron_traces
    
    return traces  # Return the list of traces
