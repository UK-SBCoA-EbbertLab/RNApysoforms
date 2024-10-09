import plotly.graph_objects as go
import polars as pl
from typing import List, Optional, Dict
from RNApysoforms.utils import check_df
import plotly.express as px

def make_traces(
    annotation: pl.DataFrame,
    expression_matrix: pl.DataFrame,
    order_transcripts_by_expression_matrix: bool = True,
    y: str = "transcript_id",
    x_start: str = "start",
    x_end: str = "end",
    annotation_hue = None,
    expression_hue = None,
    cds: str = "CDS",
    exon: str = "exon",
    intron: str = "intron",
    expression_columns: List[str] = ["counts"],
    sample_id_column: str = "sample_id",
    annotation_fill_color="grey",
    expression_fill_color="grey",
    annotation_color_palette: List[str] = px.colors.qualitative.Plotly,
    expression_color_palette: List[str] = px.colors.qualitative.Plotly_r,
    annotation_color_map: dict = None,  # Optional color map for hues    
    expression_color_map: dict = None,  # Optional color map for hues    
    intron_line_width: float = 0.5,
    exon_line_width: float = 0.25,
    expression_line_width: float = 0.25,
    line_color: str = "black",
    show_expression_points = True,
    expression_plot_style = "boxplot",
    expression_plot_opacity = 0.7,
    transcript_plot_opacity: float = 1,
    exon_height: float = 0.3,
    cds_height: float = 0.5,
    arrow_height: float = 0.5,
    arrow_length: float = 1,
    hover_start: str = "start",
    hover_end: str = "end"
) -> List[go.Box]:
    
    ## Make sure expression columns is in list form
    if isinstance(expression_columns, str):
        expression_columns = [expression_columns]


    # Check if both dataframes are None
    if annotation is None and expression_matrix is None:
        raise ValueError("At least one of 'annotation' or 'expression_matrix' must be provided.")
    

    if annotation is not None:
        if not isinstance(annotation, pl.DataFrame):
            raise TypeError(
                f"Expected annotation to be of type pl.DataFrame, got {type(annotation)}. "
                "You can use polars_df = pl.from_pandas(pandas_df) to convert a pandas DataFrame into a Polars DataFrame."
            )
        if annotation_hue == None:
            check_df(annotation, [y, x_start, x_end, "strand", "seqnames", hover_start, hover_end])
        else:
            check_df(annotation, [y, x_start, x_end, "strand", "seqnames", hover_start, hover_end, annotation_hue])
    else:
        order_transcripts_by_expression_matrix = True
        
    if expression_matrix is not None:
        if not isinstance(expression_matrix, pl.DataFrame):
            raise TypeError(
                f"Expected expression_matrix to be of type pl.DataFrame, got {type(expression_matrix)}. "
                "You can use polars_df = pl.from_pandas(pandas_df) to convert a pandas DataFrame into a Polars DataFrame."
            )
        
        required_columns = [y, sample_id_column] + expression_columns
        if expression_hue == None:
            check_df(expression_matrix, required_columns)
        else:
            required_columns.append(expression_hue)
            check_df(expression_matrix, required_columns)
    else:
        order_transcripts_by_expression_matrix = False
        
    if ((expression_matrix is not None) and (annotation is not None)):
        if not (sorted(expression_matrix["transcript_id"].unique().to_list()) == sorted(annotation["transcript_id"].unique().to_list())):
            raise ValueError("Expression matrix and annotation must contain the same transcript ids")
        

    ## Get ordering right:
    if order_transcripts_by_expression_matrix:
        # Combine transcripts
        unique_transcripts = expression_matrix[y].unique(maintain_order=True).to_list()
        # Create y_dict mapping transcript IDs to y positions
        y_dict = {val: i for i, val in enumerate(unique_transcripts)}
        # Create the sorting expression using match_literal
        annotation = annotation.with_columns(
            pl.col(y).cast(pl.Categorical).cast(pl.Utf8).replace(
                {k: i for i, k in enumerate(unique_transcripts)},
                default=len(unique_transcripts)  # Items not in custom_order will be placed at the end
            ).alias("sort_key")
        ).sort("sort_key").drop("sort_key")

    else:
        # Combine transcripts
        unique_transcripts = annotation[y].unique(maintain_order=True).to_list()
        # Create y_dict mapping transcript IDs to y positions
        y_dict = {val: i for i, val in enumerate(unique_transcripts)}
        # Create the sorting expression using match_literal
        expression_matrix = expression_matrix.with_columns(
            pl.col(y).cast(pl.Categorical).cast(pl.Utf8).replace(
                {k: i for i, k in enumerate(unique_transcripts)},
                default=len(unique_transcripts)  # Items not in custom_order will be placed at the end
            ).alias("sort_key")
        ).sort("sort_key").drop("sort_key")
        
    # Generate a color map if not provided and 'hue' is specified
    if (annotation_color_map is None) and (annotation is not None):
        values_to_colormap = annotation[annotation_hue].unique(maintain_order=True).to_list()
        annotation_color_map = {value: color for value, color in zip(values_to_colormap, annotation_color_palette)}
    elif annotation_hue is None:
        annotation_color_map = annotation_fill_color

    # Generate a color map if not provided and 'hue' is specified
    if (expression_color_map is None) and (expression_matrix is not None):
        values_to_colormap = expression_matrix[expression_hue].unique(maintain_order=True).to_list()
        expression_color_map = {value: color for value, color in zip(values_to_colormap, expression_color_palette)}
    elif expression_hue is None:
        expression_color_map = expression_fill_color

    transcript_traces = []

    # Process annotation if provided
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
        
            # Determine the fill color and legend name based on 'hue'
            if annotation_hue is None:
                exon_and_cds_color = annotation_fill_color
                hue_name = "Exon and/or CDS"
            else:
                exon_and_cds_color = annotation_color_map.get(row[annotation_hue], annotation_fill_color)
                hue_name = row[annotation_hue]
        
            # Determine whether to display the legend entry for this hue value
            if (hue_name in displayed_hue_names) or (annotation_hue is None):
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
                    hoveron='fills+points'
                )
                exon_traces.append(trace)

                if hue_name not in displayed_hue_names:
                    displayed_hue_names.append(hue_name)

            elif row["type"] == cds:
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
                    hoveron='fills+points'
                )
                cds_traces.append(trace)
        
                if hue_name not in displayed_hue_names:
                    displayed_hue_names.append(hue_name)

            elif row["type"] == intron:
                x_intron = [row[x_start], row[x_end]]
                y_intron = [y_pos, y_pos]
        
                # Add directional arrows for introns if they are sufficiently long
                if abs(row[x_start] - row[x_end]) > size / 15:
                    arrow_x = (row[x_start] + row[x_end]) / 2
                    arrow_length_px = size / (150 / arrow_length) if arrow_length != 0 else 0  
        
                    if row["strand"] == "+":  # Positive strand (arrow pointing left)
                        x_arrow = [
                            None,
                            arrow_x, arrow_x - arrow_length_px, None,
                            arrow_x, arrow_x - arrow_length_px, None
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
        
                    x_combined = x_intron + x_arrow
                    y_combined = y_intron + y_arrow
                else:
                    x_combined = x_intron
                    y_combined = y_intron

                # Create the scatter trace for the intron with optional arrows
                trace = dict(
                    type='scatter',
                    mode='lines',
                    x=x_combined,
                    y=y_combined,
                    line=dict(color=line_color, width=intron_line_width),
                    opacity=1,
                    hovertemplate=hovertemplate_text,
                    showlegend=False
                )
                intron_traces.append(trace)

        # Combine all traces (exons, CDS, introns)
        transcript_traces.extend(exon_traces + cds_traces + intron_traces)
        transcript_traces = [transcript_traces]


    # Process expression_matrix if provided
    if expression_matrix is not None:

        show_legend = True
        expression_traces = []

        # Map transcript IDs to y positions
        expression_matrix = expression_matrix.with_columns(
            pl.col(y).replace(y_dict).alias("y_pos")
        )

        if expression_hue is not None:

            # List all unique hue values
            unique_hues = expression_matrix[expression_hue].unique(maintain_order=True).to_list()
                    
            # Iterate over expression columns
            for x in expression_columns:

                x_traces_list=[]
                # Iterate over each unique hue to create Box traces
                for hue_val in unique_hues:
                    hue_filtered_df = expression_matrix.filter(pl.col(expression_hue) == hue_val)

                    box_trace = go.Box(
                        y=hue_filtered_df["y_pos"].to_list(),
                        x=hue_filtered_df[x].to_list(),
                        text=hue_filtered_df[sample_id_column].to_list(),
                        name=str(hue_val),
                        boxpoints='all',
                        jitter=0.3,
                        pointpos=0,
                        marker=dict(color=expression_color_map[hue_val]),
                        line=dict(width=expression_line_width),
                        fillcolor=expression_color_map[hue_val],
                        boxmean=False,
                        orientation='h',
                        legendgroup=str(hue_val),
                        showlegend=show_legend
                    )
                    x_traces_list.append(box_trace)

                show_legend = False
                expression_traces.append(x_traces_list)

        else:
            unique_transcripts = expression_matrix[y].unique(maintain_order=True).to_list()
            
            # Map transcript IDs to y positions
            expression_matrix = expression_matrix.with_columns(
                pl.col(y).replace(y_dict).alias("y_pos")
            )
            
            # Iterate over expression columns
            for x in expression_columns:
                x_traces_list=[]
                for transcript in unique_transcripts:
                    transcript_df = expression_matrix.filter(pl.col(y) == transcript)
                    expression = transcript_df[x].to_list()
                    sample_id = transcript_df[sample_id_column].to_list()
                    y_pos = y_dict[transcript]
                    
                    box_trace = go.Box(
                        y=[y_pos]*len(expression),
                        x=expression,
                        text=sample_id,
                        name=transcript,
                        pointpos=0,
                        boxmean=False,
                        boxpoints='all',
                        marker_color="black",
                        line=dict(width=expression_line_width),
                        fillcolor=expression_fill_color,
                        orientation='h',
                        showlegend=show_legend
                    )
                    x_traces_list.extend(box_trace)

                show_legend = False
                expression_traces.append(x_traces_list)
    traces = []

    if annotation is not None:
        traces.extend(transcript_traces)
    if expression_matrix is not None:
        traces.extend(expression_traces)

    traces.append(y_dict)

    return traces