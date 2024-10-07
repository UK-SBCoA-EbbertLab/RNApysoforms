import plotly.graph_objects as go
from plotly.subplots import make_subplots
import polars as pl
from typing import List, Tuple, Optional
from RNApysoforms.utils import check_df  # Ensure this is correctly imported
from RNApysoforms import set_axis

def make_plot(
    transcript_structure_traces: List[go.Trace],
    annotation: pl.DataFrame,
    expression_traces: Optional[List[List[go.Trace]]] = None,
    subplot_titles: List[str] = ["Transcript Structure"],
    horizontal_spacing: float = 0.1,
    vertical_spacing: float = 0.1,
    showlegend: bool = True,
    height: int = 800,
    width: int = 1800,
    template: str = "plotly_white",
) -> go.Figure:
    """
    Create a Plotly figure with three subplots: Transcript Structure, Counts, and Relative Abundance.

    Parameters:
    - transcript_structure_traces (List[go.Trace]): List of Plotly traces for Transcript Structure.
    - counts_df (pl.DataFrame): Polars DataFrame containing counts data.
    - relative_abundance_df (pl.DataFrame, optional): Polars DataFrame containing relative abundance data.
      If None, the function will generate it by passing `counts_df` with x="relative_abundance" to `make_expression_traces_func`.
    - rescaled_annotation (pl.DataFrame, optional): Polars DataFrame containing annotation data for setting genomic axis.
    - make_expression_traces_func (Callable, optional): Function to create expression traces.
      Expected signature: make_expression_traces_func(df: pl.DataFrame, x: str) -> List[go.Trace]
      Defaults to `pt.make_transcript_expression_traces`.
    - set_axis_func (Callable, optional): Function to set axis ranges or configurations.
      Expected signature: set_axis_func(fig: go.Figure, data: pl.DataFrame, x_start: str, x_end: str, row: int, col: int) -> go.Figure
      Defaults to `pt.set_axis`.
    - subplot_titles (Tuple[str, str, str], optional): Titles for each subplot.
    - horizontal_spacing (float, optional): Horizontal spacing between subplots.
    - vertical_spacing (float, optional): Vertical spacing between subplots.
    - showlegend (bool, optional): Whether to display the legend.
    - height (int, optional): Height of the figure in pixels.
    - width (int, optional): Width of the figure in pixels.
    - template (str, optional): Plotly template for styling.
    - axis_titles (dict, optional): Dictionary to customize axis titles.
      Example:
          {
              "Transcript Structure": {"x": "Genomic Position", "y": "Genes"},
              "Counts": {"x": "Transcript ID", "y": "Counts"},
              "Relative Abundance": {"x": "Transcript ID", "y": "Relative Abundance"}
          }

    Returns:
    - go.Figure: The finalized Plotly figure with three subplots.
    """
    if expression_traces != None:
        full_trace_list = transcript_structure_traces + expression_traces
    else:
        full_trace_list = [transcript_structure_traces]

    # Initialize the subplot figure with 1 row and dynamic number of columns
    fig = make_subplots(
        rows=1,
        cols=len(full_trace_list),
        subplot_titles=subplot_titles,
        horizontal_spacing=horizontal_spacing,
        vertical_spacing=vertical_spacing
    )

    for i in range(len(full_trace_list)):

        # Batch Add Transcript Structure Traces to the First Subplot
        fig.add_traces(full_trace_list[i],
            rows=[1] * len(full_trace_list[i]),
            cols=[i+1] * len(full_trace_list[i]))


    fig.update_layout(
        row=1,
        col=1,
        yaxis=dict(
            tickvals=list(range(annotation.n_unique(subset="transcript_id"))),               # Positions of the ticks
            ticktext=annotation.select("transcript_id").unique(maintain_order=True).to_series().to_list(),  # Custom labels for the ticks
            tickfont=dict(size=10, family='DejaVu Sans', color='black')),
        xaxis=dict(showticklabels=False),
        legend=dict(traceorder="reversed"),
        hoverlabel=dict(font=dict(size=8)),
        showlegend=True,
        yaxis_title="",
        xaxis_title="") 
    
    for i in range(2, (len(expression_traces) + 2)):
        fig.update_layout(
            row=1,
            col=i,
            yaxis=dict(showticklabels=False),
            legend=dict(traceorder="reversed"),
            hoverlabel=dict(font=dict(size=8)),
            showlegend=True,
            yaxis_title="",
            xaxis_title="") 

    # Update Overall Layout
    fig.update_layout(
        title_text="Transcript Analysis",
        showlegend=showlegend,
        height=height,
        width=width,
        template=template
    )

    return fig
