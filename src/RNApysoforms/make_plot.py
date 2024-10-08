import plotly.graph_objects as go
from plotly.subplots import make_subplots
import polars as pl
from typing import List, Optional
from RNApysoforms.utils import check_df  # Ensure this is correctly imported
from RNApysoforms import set_axis

def make_plot(
    transcript_structure_traces: List[go.Trace],
    annotation: pl.DataFrame,
    expression_traces: Optional[List[List[go.Trace]]] = None,
    y: str = "transcript_id",
    subplot_titles: List[str] = ["Transcript Structure"],
    horizontal_spacing: float = 0.02,
    vertical_spacing: float = 0.02,
    showlegend: bool = True,
    height: int = 800,
    width: int = 1800,
    template: str = "plotly_white",
) -> go.Figure:
    """
    Create a Plotly figure with multiple aligned subplots.
    
    [Docstring remains unchanged for brevity]
    """
    if expression_traces is not None:
        full_trace_list = [transcript_structure_traces] + expression_traces
        # Update subplot titles if needed
        if len(subplot_titles) < len(full_trace_list):
            additional_titles = [f"Expression Metric {i}" for i in range(2, len(full_trace_list)+1)]
            subplot_titles += additional_titles
    else:
        full_trace_list = [transcript_structure_traces]

    # Initialize the subplot figure with shared y-axes
    fig = make_subplots(
        rows=1,
        cols=len(full_trace_list),
        subplot_titles=subplot_titles,
        horizontal_spacing=horizontal_spacing,
        vertical_spacing=vertical_spacing,
        shared_yaxes=True  # Share y-axes across all subplots
    )


    # Create a consistent y-axis mapping
    unique_transcripts = annotation[y].unique(maintain_order=True).to_list()
    y_dict = {val: i for i, val in enumerate(unique_transcripts)}

    # Add all traces to their respective subplots
    for i, subplot_traces in enumerate(full_trace_list, start=1):
        fig.add_traces(
            subplot_traces,
            rows=[1] * len(subplot_traces),
            cols=[i] * len(subplot_traces)
        )

    # Customize the shared y-axis (only for the first subplot)
    fig.update_yaxes(
        tickvals=list(y_dict.values()),
        ticktext=list(y_dict.keys()),
        tickfont=dict(size=10, family='DejaVu Sans', color='black'),
        title="",  # Optional
        row=1,
        col=1,
        range=[-0.8, (len(unique_transcripts) - 0.2)]
    )

    # Customize x-axes
    # First subplot (Transcript Structure)
    fig.update_xaxes(
        showticklabels=False,
        title="",
        row=1,
        col=1
    )

    # Additional subplots (Expression Traces)
    if expression_traces:
        for i in range(2, len(full_trace_list) + 1):
            fig.update_xaxes(
                showticklabels=True,  # Show x-axis labels for expression plots
                title=f"",  # Customize as needed
                row=1,
                col=i
            )
            # Hide y-axis labels for additional subplots
            fig.update_yaxes(
                tickvals=list(y_dict.values()),
                ticktext=list(y_dict.keys()),
                showticklabels=False,
                row=1,
                col=i,
                range=[-0.8, (len(unique_transcripts) - 0.2)]
            )


    # Update overall layout
    fig.update_layout(
        title_text="",
        showlegend=showlegend,
        height=height,
        width=width,
        template=template,
        legend=dict(traceorder="reversed"),
        hoverlabel=dict(font=dict(size=12)),
        margin=dict(l=100, r=50, t=100, b=50),  # Adjust margins as needed
        boxmode='group'
    )

    return fig