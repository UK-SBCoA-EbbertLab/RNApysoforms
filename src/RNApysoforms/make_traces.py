import plotly.graph_objects as go
import polars as pl
from typing import List, Optional, Dict, Union
from RNApysoforms.utils import check_df
import plotly.express as px
import warnings

def make_traces(
    annotation: Optional[pl.DataFrame] = None,
    expression_matrix: Optional[pl.DataFrame] = None,
    order_transcripts_by_expression_matrix: bool = True,
    y: str = "transcript_id",
    x_start: str = "start",
    x_end: str = "end",
    annotation_hue: Optional[str] = None,
    expression_hue: Optional[str] = None,
    cds: str = "CDS",
    exon: str = "exon",
    intron: str = "intron",
    expression_columns: Union[str, List[str]] = ["counts"],
    sample_id_column: str = "sample_id",
    annotation_fill_color: str = "grey",
    expression_fill_color: str = "grey",
    annotation_color_palette: List[str] = px.colors.qualitative.Plotly,
    expression_color_palette: List[str] = px.colors.qualitative.Plotly_r,
    annotation_color_map: Optional[Dict[str, str]] = None,  # Optional color map for hues
    expression_color_map: Optional[Dict[str, str]] = None,  # Optional color map for hues
    intron_line_width: float = 0.5,
    exon_line_width: float = 0.25,
    expression_line_width: float = 0.5,
    line_color: str = "black",
    expression_plot_style: str = "boxplot",
    spanmode: str = "hard",
    marker_color: str = "black",
    marker_opacity: float = 1,
    marker_size: int = 5,
    marker_jitter: float = 0.3,
    expression_plot_opacity: float = 1,
    transcript_plot_opacity: float = 1,
    exon_height: float = 0.3,
    cds_height: float = 0.5,
    arrow_height: float = 0.5,  # Note: Not currently used in the function
    arrow_length: float = 1,    # Note: Not currently used in the function
    arrow_size: float = 10,
    hover_start: str = "start",
    hover_end: str = "end",
    show_box_mean: bool = True,
    box_points: Union[str, bool] = "all",
    expression_plot_legend_title: str = "<b><u>Expression Plot Hue<u><b>",
    transcript_plot_legend_title: str = "<b><u>Transcript Structure Hue<u><b>",
) -> List[Union[go.Box, go.Violin, dict, Dict[str, int]]]:
    """
    Generates Plotly traces for visualizing transcript structures and expression data.

    This function processes genomic annotation and expression data to create Plotly traces suitable for plotting transcript
    structures alongside expression data. It supports customization of plot aesthetics, including colors, line widths,
    plot styles, and annotations. The function returns a list of traces that can be directly used in Plotly figures,
    including traces for exons, introns, CDS regions, and expression data.

    **Required Columns in `annotation` DataFrame:**
    - `y` (default `"transcript_id"`): Identifier for each transcript.
    - `x_start` (default `"start"`): Start position of the feature.
    - `x_end` (default `"end"`): End position of the feature.
    - `"strand"`: Strand information (`"+"` or `"-"`).
    - `"seqnames"`: Chromosome or sequence name.
    - `hover_start` (default `"start"`): Start position for hover information.
    - `hover_end` (default `"end"`): End position for hover information.
    - If `annotation_hue` is specified, it must also be a column in `annotation`.

    **Required Columns in `expression_matrix` DataFrame:**
    - `y` (default `"transcript_id"`): Identifier for each transcript.
    - `sample_id_column` (default `"sample_id"`): Identifier for each sample.
    - `expression_columns` (default `["counts"]`): Column name or list of column names containing expression values
      that you want to plot in the order you want to plot them.
    - If `expression_hue` is specified, it must also be a column in `expression_matrix`.

    Parameters
    ----------
    annotation : pl.DataFrame, optional
        A Polars DataFrame containing genomic annotation data for transcripts. Includes exons, introns, and CDS features.
        If provided, the function will generate traces for transcript structures.
    expression_matrix : pl.DataFrame, optional
        A Polars DataFrame containing expression data. If provided, the function will generate traces for expression plots.
    order_transcripts_by_expression_matrix : bool, optional
        If True, orders transcripts based on their order in the expression matrix. If False, orders by annotation DataFrame.
        Default is True.
    y : str, optional
        Column name in both `annotation` and `expression_matrix` representing transcript identifiers. Default is "transcript_id".
    x_start : str, optional
        Column name in `annotation` representing the start position of features. Default is "start".
    x_end : str, optional
        Column name in `annotation` representing the end position of features. Default is "end".
    annotation_hue : str, optional
        Column name in `annotation` used to color-code transcript features based on categories. Default is None.
    expression_hue : str, optional
        Column name in `expression_matrix` used to color-code expression data based on categories. Default is None.
    cds : str, optional
        Value in the `type` column of `annotation` representing CDS features. Default is "CDS".
    exon : str, optional
        Value in the `type` column of `annotation` representing exon features. Default is "exon".
    intron : str, optional
        Value in the `type` column of `annotation` representing intron features. Default is "intron".
    expression_columns : Union[str, List[str]], optional
        Column name or list of column names in `expression_matrix` containing expression values. If a string is provided,
        it is converted to a list containing that string. Default is ["counts"].
    sample_id_column : str, optional
        Column name in `expression_matrix` representing sample identifiers. Default is "sample_id".
    annotation_fill_color : str, optional
        Default fill color for transcript features if `annotation_hue` is not specified. Default is "grey".
    expression_fill_color : str, optional
        Default fill color for expression plots if `expression_hue` is not specified. Default is "grey".
    annotation_color_palette : List[str], optional
        List of colors to use for different categories in `annotation_hue`. Default uses Plotly qualitative palette.
    expression_color_palette : List[str], optional
        List of colors to use for different categories in `expression_hue`. Default uses reversed Plotly qualitative palette.
    annotation_color_map : dict, optional
        Mapping from categories in `annotation_hue` to colors. If None, colors are assigned from `annotation_color_palette`.
    expression_color_map : dict, optional
        Mapping from categories in `expression_hue` to colors. If None, colors are assigned from `expression_color_palette`.
    intron_line_width : float, optional
        Line width for intron traces. Default is 0.5.
    exon_line_width : float, optional
        Line width for exon traces. Default is 0.25.
    expression_line_width : float, optional
        Line width for expression plot traces. Default is 0.5.
    line_color : str, optional
        Color for the lines outlining transcript features. Default is "black".
    expression_plot_style : str, optional
        Style of the expression plot. Options are "boxplot" or "violin". Default is "boxplot".
    spanmode : str, optional
        For violin plots, defines how the width of the violin spans the data. Options are "hard" or "soft". Default is "hard".
    marker_color : str, optional
        Color of the markers in expression plots. Default is "black".
    marker_opacity : float, optional
        Opacity of the markers in expression plots. Default is 1.
    marker_size : int, optional
        Size of the markers in expression plots. Default is 5.
    marker_jitter : float, optional
        Amount of jitter (spread) applied to the markers in expression plots. Default is 0.3.
    expression_plot_opacity : float, optional
        Opacity of the expression plot traces. Default is 1.
    transcript_plot_opacity : float, optional
        Opacity of the transcript structure traces. Default is 1.
    exon_height : float, optional
        Height of exon rectangles in the plot. Default is 0.3.
    cds_height : float, optional
        Height of CDS rectangles in the plot. Default is 0.5.
    arrow_height : float, optional
        Height of the arrow markers for introns. **Note:** Currently not used in the function. Default is 0.5.
    arrow_length : float, optional
        Length of the arrow markers for introns. **Note:** Currently not used in the function. Default is 1.
    arrow_size : float, optional
        Size of the arrow markers for introns. Default is 10.
    hover_start : str, optional
        Column name in `annotation` for the start position displayed in hover information. Default is "start".
    hover_end : str, optional
        Column name in `annotation` for the end position displayed in hover information. Default is "end".
    show_box_mean : bool, optional
        If True, shows the mean in box and violin plots. Default is True.
    box_points : Union[str, bool], optional
        Controls the display of points in box and violin plots. Options include "all", "outliers", "suspectedoutliers", or False. Default is "all".
    expression_plot_legend_title : str, optional
        Title for the legend of the expression plot. Default is "<b><u>Expression Plot Hue<u><b>".
    transcript_plot_legend_title : str, optional
        Title for the legend of the transcript structure plot. Default is "<b><u>Transcript Structure Hue<u><b>".

    Returns
    -------
    List[Union[go.Box, go.Violin, dict, Dict[str, int]]]
        A list containing the generated Plotly traces and a mapping of transcript identifiers to y-axis positions.
        The list includes:
        - Transcript feature traces (exons, CDS, introns) as dictionaries.
        - Expression plot traces as `go.Box` or `go.Violin` objects.
        - The `y_dict` mapping, which maps transcript identifiers to their corresponding y-axis positions.

    Raises
    ------
    ValueError
        If neither `annotation` nor `expression_matrix` is provided.
    TypeError
        If `annotation` or `expression_matrix` are not Polars DataFrames.
    ValueError
        If required columns are missing in `annotation` or `expression_matrix`.
    ValueError
        If there are no matching transcripts between `annotation` and `expression_matrix`.
    ValueError
        If an invalid `expression_plot_style` is provided.

    Examples
    --------
    Generate traces for plotting transcript structures and expression data:

    >>> import polars as pl
    >>> from RNApysoforms.plotting import make_traces
    >>> # Prepare annotation DataFrame
    >>> annotation_df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
    ...     "start": [100, 200, 150, 250],
    ...     "end": [150, 250, 200, 300],
    ...     "type": ["exon", "CDS", "exon", "CDS"],
    ...     "strand": ["+", "+", "-", "-"],
    ...     "seqnames": ["chr1", "chr1", "chr2", "chr2"]
    ... })
    >>> # Prepare expression matrix
    >>> expression_df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
    ...     "sample_id": ["sample1", "sample2", "sample1", "sample2"],
    ...     "counts": [100, 200, 150, 250]
    ... })
    >>> # Generate traces
    >>> traces = make_traces(annotation=annotation_df, expression_matrix=expression_df)
    >>> # Use traces to create a Plotly figure (not shown here)

    Notes
    -----
    - The function ensures that both `annotation` and `expression_matrix` contain common transcripts.
      It filters out any transcripts not present in both DataFrames.
    - Warnings are issued if transcripts are present in one DataFrame but missing in the other.
    - Traces are generated for exons, CDS, and introns with customizable aesthetics.
    - Expression data can be visualized using box plots or violin plots, with options for coloring by categories.
    - The `y_dict` mapping is used to align transcripts across different plots by assigning consistent y-axis positions.
    - The function handles strand direction when plotting intron arrows.
    - Custom legends and hover information can be configured via parameters.
    - **Note:** The parameters `arrow_height` and `arrow_length` are currently not utilized in the function.

    """

    # Ensure that expression_columns is a list
    if isinstance(expression_columns, str):
        expression_columns = [expression_columns]

    # Check if both dataframes are None
    if annotation is None and expression_matrix is None:
        raise ValueError("At least one of 'annotation' or 'expression_matrix' must be provided.")

    # Validate and process the 'annotation' DataFrame
    if annotation is not None:
        # Check if 'annotation' is a Polars DataFrame
        if not isinstance(annotation, pl.DataFrame):
            raise TypeError(
                f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}. "
                "You can use pl.from_pandas(pandas_df) to convert a pandas DataFrame into a Polars DataFrame."
            )
        # Ensure required columns are present in 'annotation'
        required_columns = [y, x_start, x_end, "strand", "seqnames", hover_start, hover_end]
        if annotation_hue is not None:
            required_columns.append(annotation_hue)
        check_df(annotation, required_columns)
    else:
        # If 'annotation' is None, transcripts will be ordered by 'expression_matrix'
        order_transcripts_by_expression_matrix = True

    # Validate and process the 'expression_matrix' DataFrame
    if expression_matrix is not None:
        # Check if 'expression_matrix' is a Polars DataFrame
        if not isinstance(expression_matrix, pl.DataFrame):
            raise TypeError(
                f"Expected 'expression_matrix' to be of type pl.DataFrame, got {type(expression_matrix)}. "
                "You can use pl.from_pandas(pandas_df) to convert a pandas DataFrame into a Polars DataFrame."
            )
        # Ensure required columns are present in 'expression_matrix'
        required_columns = [y, sample_id_column] + expression_columns
        if expression_hue is not None:
            required_columns.append(expression_hue)
        check_df(expression_matrix, required_columns)
    else:
        # If 'expression_matrix' is None, transcripts will be ordered by 'annotation'
        order_transcripts_by_expression_matrix = False

    # If both 'annotation' and 'expression_matrix' are provided, ensure they have common transcripts
    if annotation is not None and expression_matrix is not None:
        # Filter 'expression_matrix' to include only transcripts present in 'annotation'
        filtered_expression_matrix = expression_matrix.filter(
            pl.col(y).is_in(annotation[y])
        )
        # Filter 'annotation' to include only transcripts present in 'expression_matrix'
        filtered_annotation = annotation.filter(
            pl.col(y).is_in(expression_matrix[y])
        )
        # Check if filtered data is empty and raise an error if true
        if filtered_expression_matrix.is_empty() or filtered_annotation.is_empty():
            raise ValueError(
                f"No matching '{y}' entries between annotation and expression matrix."
            )
        # Identify discrepancies between transcripts in 'annotation' and 'expression_matrix'
        annotation_transcripts = set(annotation[y].unique())
        expression_transcripts = set(expression_matrix[y].unique())
        missing_in_expression = annotation_transcripts - expression_transcripts
        missing_in_annotation = expression_transcripts - annotation_transcripts

        # Warn about transcripts missing in the expression matrix
        if missing_in_expression:
            warnings.warn(
                f"{len(missing_in_expression)} transcript(s) are present in the annotation but missing in the expression matrix. "
                f"Missing transcripts: {', '.join(sorted(missing_in_expression))}. "
                "Only transcripts present in both will be used for making traces."
            )
        # Warn about transcripts missing in the annotation
        if missing_in_annotation:
            warnings.warn(
                f"{len(missing_in_annotation)} transcript(s) are present in the expression matrix but missing in the annotation. "
                f"Missing transcripts: {', '.join(sorted(missing_in_annotation))}. "
                "Only transcripts present in both will be used for making traces."
            )
        # Keep only common transcripts in both DataFrames
        common_transcripts = annotation_transcripts & expression_transcripts
        annotation = annotation.filter(
            pl.col(y).is_in(common_transcripts)
        )
        expression_matrix = expression_matrix.filter(
            pl.col(y).is_in(common_transcripts)
        )

    # Determine the ordering of transcripts
    if order_transcripts_by_expression_matrix:
        # Order transcripts based on 'expression_matrix'
        unique_transcripts = expression_matrix[y].unique(maintain_order=True).to_list()
        # Create y_dict mapping transcript IDs to y positions
        y_dict = {val: i for i, val in enumerate(unique_transcripts)}
        if annotation is not None:
            # Sort 'annotation' DataFrame based on transcript order
            annotation = annotation.with_columns(
                pl.col(y).cast(pl.Categorical).cast(pl.Utf8).replace(
                    {k: i for i, k in enumerate(unique_transcripts)},
                    default=len(unique_transcripts)  # Items not in custom_order will be placed at the end
                ).alias("sort_key")
            ).sort("sort_key").drop("sort_key")
    else:
        # Order transcripts based on 'annotation'
        unique_transcripts = annotation[y].unique(maintain_order=True).to_list()
        # Create y_dict mapping transcript IDs to y positions
        y_dict = {val: i for i, val in enumerate(unique_transcripts)}
        if expression_matrix is not None:
            # Sort 'expression_matrix' DataFrame based on transcript order
            expression_matrix = expression_matrix.with_columns(
                pl.col(y).cast(pl.Categorical).cast(pl.Utf8).replace(
                    {k: i for i, k in enumerate(unique_transcripts)},
                    default=len(unique_transcripts)  # Items not in custom_order will be placed at the end
                ).alias("sort_key")
            ).sort("sort_key").drop("sort_key")

    # Generate color maps if not provided and 'hue' is specified
    if annotation_color_map is None and annotation is not None and annotation_hue is not None:
        values_to_colormap = annotation[annotation_hue].unique(maintain_order=True).to_list()
        annotation_color_map = {value: color for value, color in zip(values_to_colormap, annotation_color_palette)}
    elif annotation_hue is None and annotation_color_map is None:
        # Use default fill color if no hue or color map is specified
        annotation_color_map = annotation_fill_color

    if expression_color_map is None and expression_matrix is not None and expression_hue is not None:
        values_to_colormap = expression_matrix[expression_hue].unique(maintain_order=True).to_list()
        expression_color_map = {value: color for value, color in zip(values_to_colormap, expression_color_palette)}
    elif expression_hue is None and expression_color_map is None:
        # Use default fill color if no hue or color map is specified
        expression_color_map = expression_fill_color

    transcript_traces = []

    # Process 'annotation' to create transcript structure traces
    if annotation is not None:
        # Initialize lists to store traces for different feature types
        cds_traces = []      # Stores traces for CDS (coding sequences)
        intron_traces = []   # Stores traces for introns
        exon_traces = []     # Stores traces for exons

        # Calculate the global maximum and minimum x-values (positions)
        global_max = max(
            annotation.select(pl.col(x_start).max()).item(),
            annotation.select(pl.col(x_end).max()).item()
        )
        global_min = min(
            annotation.select(pl.col(x_start).min()).item(),
            annotation.select(pl.col(x_end).min()).item()
        )
        # Calculate the total size of the x-axis range
        size = int(abs(global_max - global_min))

        # Create a list to keep track of hue values already displayed in the legend
        displayed_hue_names = []

        # Iterate over each row in the DataFrame to create traces for exons, CDS, and introns
        for row in annotation.iter_rows(named=True):
            y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript

            # Determine the fill color and legend name based on 'annotation_hue'
            if annotation_hue is None:
                exon_and_cds_color = annotation_fill_color
                hue_name = "Exon and/or CDS"
            else:
                exon_and_cds_color = annotation_color_map.get(row[annotation_hue], annotation_fill_color)
                hue_name = row[annotation_hue]

            # Determine whether to display the legend entry for this hue value
            if hue_name in displayed_hue_names:
                display_legend = False
            else:
                display_legend = True

            # Define hover template with feature type, number, start, and end positions for each row
            feature_size = abs(row[hover_end] - row[hover_start])
            hovertemplate_text = (
                f"<b>{y}:</b> {row[y]}<br>"
                f"<b>Feature Type:</b> {row['type']}<br>"
                f"<b>Feature Number:</b> {row.get('exon_number', 'N/A')}<br>"
                f"<b>Chromosome:</b> {row['seqnames']}<br>"
                f"<b>Start:</b> {row[hover_start]}<br>"
                f"<b>End:</b> {row[hover_end]}<br>"
                f"<b>Size:</b> {feature_size}<br>"
                "<extra></extra>"
            )

            # Create trace based on the feature type
            if row["type"] == exon:
                # Define coordinates for the exon rectangle
                x0 = row[x_start]
                x1 = row[x_end]
                y0 = y_pos - exon_height / 2
                y1 = y_pos + exon_height / 2

                # Create the scatter trace for the exon
                trace = dict(
                    type='scatter',
                    mode='lines+markers',
                    x=[x0, x1, x1, x0, x0],
                    y=[y0, y0, y1, y1, y0],
                    fill='toself',
                    fillcolor=exon_and_cds_color,
                    line=dict(color=line_color, width=exon_line_width),
                    marker=dict(opacity=0),
                    opacity=transcript_plot_opacity,
                    name=hue_name,
                    legendgroup=hue_name,
                    showlegend=display_legend,
                    hovertemplate=hovertemplate_text,
                    hoverlabel=dict(namelength=-1),
                    hoveron='fills+points',
                    legendgrouptitle_text=transcript_plot_legend_title
                )
                exon_traces.append(trace)
                # Prevent duplicate legend entries
                if hue_name not in displayed_hue_names:
                    displayed_hue_names.append(hue_name)
                transcript_plot_legend_title = ""  # Reset legend title after first use

            elif row["type"] == cds:
                # Define coordinates for the CDS rectangle
                x0 = row[x_start]
                x1 = row[x_end]
                y0 = y_pos - cds_height / 2
                y1 = y_pos + cds_height / 2

                # Create the scatter trace for the CDS
                trace = dict(
                    type='scatter',
                    mode='lines+markers',
                    x=[x0, x1, x1, x0, x0],
                    y=[y0, y0, y1, y1, y0],
                    fill='toself',
                    fillcolor=exon_and_cds_color,
                    line=dict(color=line_color, width=exon_line_width),
                    marker=dict(opacity=0),
                    opacity=transcript_plot_opacity,
                    name=hue_name,
                    legendgroup=hue_name,
                    showlegend=display_legend,
                    hovertemplate=hovertemplate_text,
                    hoverlabel=dict(namelength=-1),
                    hoveron='fills+points',
                    legendgrouptitle_text=transcript_plot_legend_title
                )
                cds_traces.append(trace)
                # Prevent duplicate legend entries
                if hue_name not in displayed_hue_names:
                    displayed_hue_names.append(hue_name)
                transcript_plot_legend_title = ""  # Reset legend title after first use

            elif row["type"] == intron:
                # Define coordinates for the intron line
                x_intron = [row[x_start], row[x_end]]
                y_intron = [y_pos, y_pos]

                # Add an arrow marker if the intron is sufficiently long
                if abs(row[x_start] - row[x_end]) > size / 15:
                    if row["strand"] == "-":
                        # Arrow pointing left, placed before the intron start
                        marker_symbol = 'arrow-left'
                        arrow_x = ((row[x_start] + row[x_end]) / 2) - abs((row[x_end] - row[x_start]) / 7)
                    elif row["strand"] == "+":
                        # Arrow pointing right, placed after the intron start
                        marker_symbol = 'arrow-right'
                        arrow_x = ((row[x_start] + row[x_end]) / 2) + abs((row[x_end] - row[x_start]) / 7)
                    arrow_y = y_pos

                    # Create the scatter trace for the arrow marker
                    trace_arrow = dict(
                        type='scatter',
                        mode='markers',
                        x=[arrow_x],
                        y=[arrow_y],
                        marker=dict(symbol=marker_symbol, size=arrow_size, color=line_color),
                        opacity=1,
                        hoverinfo='skip',  # Skip hover info for the arrow
                        showlegend=False
                    )
                    intron_traces.append(trace_arrow)

                # Create the scatter trace for the intron line
                trace_intron = dict(
                    type='scatter',
                    mode='lines',
                    x=x_intron,
                    y=y_intron,
                    line=dict(color=line_color, width=intron_line_width),
                    opacity=1,
                    hovertemplate=hovertemplate_text,
                    showlegend=False
                )
                intron_traces.append(trace_intron)

        # Combine all traces (exons, CDS, introns)
        transcript_traces.extend(exon_traces + cds_traces + intron_traces)
        transcript_traces = [transcript_traces]  # Wrap in a list to maintain consistency

    # Process 'expression_matrix' to create expression plot traces
    if expression_matrix is not None:
        show_legend = True  # Control legend display
        expression_traces = []

        # Map transcript IDs to y positions
        expression_matrix = expression_matrix.with_columns(
            pl.col(y).replace(y_dict).alias("y_pos")
        )

        if expression_hue is not None:
            # List all unique hue values
            unique_hues = expression_matrix[expression_hue].unique(maintain_order=True).to_list()
            # Create offset mapping for grouped plots
            offset_dict = {val: i for i, val in enumerate(unique_hues)}

            # Iterate over expression columns
            for x in expression_columns:
                x_traces_list = []
                # Iterate over each unique hue to create traces
                for hue_val in unique_hues:
                    hue_filtered_df = expression_matrix.filter(pl.col(expression_hue) == hue_val)
                    # Create the appropriate plot based on 'expression_plot_style'
                    if expression_plot_style == "boxplot":
                        box_trace = go.Box(
                            y=hue_filtered_df["y_pos"].to_list(),
                            x=hue_filtered_df[x].to_list(),
                            text=hue_filtered_df[sample_id_column].to_list(),
                            name=str(hue_val),
                            boxpoints=box_points,
                            jitter=marker_jitter,
                            pointpos=0,
                            line=dict(width=expression_line_width),
                            fillcolor=expression_color_map[hue_val],
                            boxmean=show_box_mean,
                            orientation='h',
                            legendgroup=str(hue_val),
                            showlegend=show_legend,
                            offsetgroup=offset_dict[hue_val],
                            opacity=expression_plot_opacity,
                            marker=dict(opacity=marker_opacity, size=marker_size, color=marker_color),
                            legendgrouptitle_text=expression_plot_legend_title,
                        )
                    elif expression_plot_style == "violin":
                        box_trace = go.Violin(
                            y=hue_filtered_df["y_pos"].to_list(),
                            x=hue_filtered_df[x].to_list(),
                            text=hue_filtered_df[sample_id_column].to_list(),
                            name=str(hue_val),
                            points=box_points,
                            jitter=marker_jitter,
                            pointpos=0,
                            line=dict(width=expression_line_width),
                            fillcolor=expression_color_map[hue_val],
                            meanline_visible=show_box_mean,
                            orientation='h',
                            legendgroup=str(hue_val),
                            showlegend=show_legend,
                            offsetgroup=offset_dict[hue_val],
                            opacity=expression_plot_opacity,
                            marker=dict(opacity=marker_opacity, size=marker_size, color=marker_color),
                            legendgrouptitle_text=expression_plot_legend_title,
                            spanmode=spanmode
                        )
                    else:
                        raise ValueError(f"Invalid expression_plot_style: {expression_plot_style}")
                    x_traces_list.append(box_trace)
                    expression_plot_legend_title = ""  # Reset legend title after first use
                show_legend = False  # Only show legend once
                expression_traces.append(x_traces_list)

        else:
            # No 'expression_hue' specified
            unique_transcripts = expression_matrix[y].unique(maintain_order=True).to_list()
            # Map transcript IDs to y positions
            expression_matrix = expression_matrix.with_columns(
                pl.col(y).replace(y_dict).alias("y_pos")
            )

            # Iterate over expression columns
            for x in expression_columns:
                x_traces_list = []
                for transcript in unique_transcripts:
                    transcript_df = expression_matrix.filter(pl.col(y) == transcript)
                    expression = transcript_df[x].to_list()
                    sample_id = transcript_df[sample_id_column].to_list()
                    y_pos = y_dict[transcript]

                    # Create the appropriate plot based on 'expression_plot_style'
                    if expression_plot_style == "boxplot":
                        box_trace = go.Box(
                            y=[y_pos] * len(expression),
                            x=expression,
                            text=sample_id,
                            name="Box Plots",
                            pointpos=0,
                            offsetgroup=0,
                            boxmean=show_box_mean,
                            jitter=marker_jitter,
                            boxpoints=box_points,
                            line=dict(width=expression_line_width),
                            fillcolor=expression_fill_color,
                            orientation='h',
                            showlegend=show_legend,
                            opacity=expression_plot_opacity,
                            marker=dict(opacity=marker_opacity, size=marker_size, color=marker_color),
                            legendgrouptitle_text=expression_plot_legend_title,
                            legendgroup="expression"
                        )
                    elif expression_plot_style == "violin":
                        box_trace = go.Violin(
                            y=[y_pos] * len(expression),
                            x=expression,
                            text=sample_id,
                            name="Violin Plots",
                            pointpos=0,
                            offsetgroup=0,
                            meanline_visible=show_box_mean,
                            jitter=marker_jitter,
                            points=box_points,
                            line=dict(width=expression_line_width),
                            fillcolor=expression_fill_color,
                            orientation='h',
                            showlegend=show_legend,
                            opacity=expression_plot_opacity,
                            marker=dict(opacity=marker_opacity, size=marker_size, color=marker_color),
                            legendgrouptitle_text=expression_plot_legend_title,
                            legendgroup="expression",
                            spanmode=spanmode
                        )
                    else:
                        raise ValueError(f"Invalid expression_plot_style: {expression_plot_style}")
                    x_traces_list.append(box_trace)
                    expression_plot_legend_title = ""
                    show_legend = False  # Only show legend once
                expression_traces.append(x_traces_list)

    # Combine transcript and expression traces
    traces = []
    if annotation is not None:
        traces.extend(transcript_traces)
    if expression_matrix is not None:
        traces.extend(expression_traces)
    # Append the y_dict mapping to the traces list
    traces.append(y_dict)

    return traces  # Return the list of traces and y-axis mapping
