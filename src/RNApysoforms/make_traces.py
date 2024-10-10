import plotly.graph_objects as go
import polars as pl
from typing import List, Optional, Dict
from RNApysoforms.utils import check_df
import plotly.express as px
import warnings

def make_traces(
    annotation: pl.DataFrame = None,
    expression_matrix: pl.DataFrame = None,
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
    expression_line_width: float = 0.5,
    line_color: str = "black",
    expression_plot_style = "boxplot",
    spanmode="hard",
    marker_color = "black",
    marker_opacity = 1,
    marker_size = 5,
    marker_jitter = 0.3,
    expression_plot_opacity = 1,
    transcript_plot_opacity: float = 1,
    exon_height: float = 0.3,
    cds_height: float = 0.5,
    arrow_height: float = 0.5,
    arrow_length: float = 1,
    arrow_size: float = 10,
    hover_start: str = "start",
    hover_end: str = "end",
    show_box_mean=True,
    box_points="all",
    expression_plot_legend_title="<b><u>Expression Plot Hue<u><b>",
    transcript_plot_legend_title="<b><u>Transcript Structure Hue<u><b>",
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
        
        # Filter expression_matrix to include only transcripts present in filtered_annotation
        filtered_expression_matrix = expression_matrix.filter(
            pl.col(y).is_in(annotation[y])
        )

        # Filter expression_matrix to include only transcripts present in filtered_annotation
        filtered_annotation = annotation.filter(
            pl.col(y).is_in(expression_matrix[y])
        )

        # If filtered expression matrix is empty, raise an error
        if (filtered_expression_matrix.is_empty()) or (filtered_annotation.is_empty()):
            raise ValueError(
                f"No matching '{y}' entries between annotation and expression matrix."
            )

        # Check for discrepancies between transcripts in annotation and expression_matrix
        annotation_transcripts = set(filtered_annotation[y].unique())
        expression_transcripts = set(filtered_expression_matrix[y].unique())

        # Transcripts in annotation but not in expression matrix
        missing_in_expression = annotation_transcripts - expression_transcripts

        # Transcripts in expression matrix but not in annotation
        missing_in_annotation = expression_transcripts - annotation_transcripts

        # Warn about transcripts missing in the expression matrix
        if missing_in_expression:
            warnings.warn(
                f"{len(missing_in_expression)} transcript(s) are present in the annotation but missing in the expression matrix. "
                f"Missing transcripts: {', '.join(sorted(missing_in_expression))}. "
                "Only transcripts present in both will be used for making traces"
            )

        # Warn about transcripts missing in the annotation
        if missing_in_annotation:
            warnings.warn(
                f"{len(missing_in_annotation)} transcript(s) are present in the expression matrix but missing in the annotation. "
                f"Missing transcripts: {', '.join(sorted(missing_in_annotation))}. "
                "Only transcripts present in both will be used for making traces"
            )

        # Ensure both filtered_annotation and filtered_expression_matrix contain only common transcripts
        common_transcripts = annotation_transcripts & expression_transcripts
        annotation = annotation.filter(
            pl.col(y).is_in(common_transcripts)
        )
        expression_matrix = expression_matrix.filter(
            pl.col(y).is_in(common_transcripts)
        )
        
        

    ## Get ordering right:
    if order_transcripts_by_expression_matrix:
        # Combine transcripts
        unique_transcripts = expression_matrix[y].unique(maintain_order=True).to_list()
        # Create y_dict mapping transcript IDs to y positions
        y_dict = {val: i for i, val in enumerate(unique_transcripts)}
        
        if annotation is not None:
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
        if expression_matrix is not None:
            # Create the sorting expression using match_literal
            expression_matrix = expression_matrix.with_columns(
                pl.col(y).cast(pl.Categorical).cast(pl.Utf8).replace(
                    {k: i for i, k in enumerate(unique_transcripts)},
                    default=len(unique_transcripts)  # Items not in custom_order will be placed at the end
                ).alias("sort_key")
            ).sort("sort_key").drop("sort_key")
            
    # Generate a color map if not provided and 'hue' is specified
    if (annotation_color_map is None) and (annotation is not None) and (annotation_hue is not None):
        values_to_colormap = annotation[annotation_hue].unique(maintain_order=True).to_list()
        annotation_color_map = {value: color for value, color in zip(values_to_colormap, annotation_color_palette)}
    elif ((annotation_hue is None) and (annotation_color_map is None)):
        annotation_color_map = annotation_fill_color

    # Generate a color map if not provided and 'hue' is specified
    if (expression_color_map is None) and (expression_matrix is not None) and (expression_hue is not None):
        values_to_colormap = expression_matrix[expression_hue].unique(maintain_order=True).to_list()
        expression_color_map = {value: color for value, color in zip(values_to_colormap, expression_color_palette)}
    elif ((expression_hue is None) and (expression_color_map is None)):
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
            if (hue_name in displayed_hue_names):
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
                    hoveron='fills+points',
                    legendgrouptitle_text=transcript_plot_legend_title
                )
                transcript_plot_legend_title=""
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
                    hoveron='fills+points',
                    legendgrouptitle_text=transcript_plot_legend_title
                )
                
                cds_traces.append(trace)
                transcript_plot_legend_title=""

                if hue_name not in displayed_hue_names:
                    displayed_hue_names.append(hue_name)

            elif row["type"] == intron:

                
                x_intron = [row[x_start], row[x_end]]
                y_intron = [y_pos, y_pos]

                # Add an arrow marker before the intron if it's sufficiently long
                if abs(row[x_start] - row[x_end]) > size / 15:
                    
                    if row["strand"] == "-":
                        # Arrow pointing left, placed before the intron start
                        marker_symbol = 'arrow-left'
                        arrow_x = ((row[x_start] + row[x_end])/2) - abs((row[x_end] - row[x_start])/7)
                    elif row["strand"] == "+":
                        # Arrow pointing right, placed before the intron end
                        marker_symbol = 'arrow-right'
                        arrow_x = ((row[x_start] + row[x_end])/2) + abs((row[x_end] - row[x_start])/7)

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
            # Create y_dict mapping transcript IDs to y positions
            offset_dict = {val: i for i, val in enumerate(unique_hues)} 
                    
            # Iterate over expression columns
            for x in expression_columns:

                x_traces_list=[]
                # Iterate over each unique hue to create Box traces
                for hue_val in unique_hues:
                    hue_filtered_df = expression_matrix.filter(pl.col(expression_hue) == hue_val)

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
                        print("error")

                    x_traces_list.append(box_trace)
                    expression_plot_legend_title = ""

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
                    
                    if expression_plot_style == "boxplot":
                        box_trace = go.Box(
                            y=[y_pos]*len(expression),
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
                            legendgroup="expression")
                    
                    elif expression_plot_style == "violin":
                        box_trace = go.Violin(
                            y=[y_pos]*len(expression),
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
                            spanmode=spanmode)
                    else:
                        print("error")    
                    
                    x_traces_list.append(box_trace)
                    expression_plot_legend_title=""
                    show_legend = False

                expression_traces.append(x_traces_list)
    traces = []

    if annotation is not None:
        traces.extend(transcript_traces)
    if expression_matrix is not None:
        traces.extend(expression_traces)

    traces.append(y_dict)

    return traces