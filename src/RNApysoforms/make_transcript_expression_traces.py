import plotly.graph_objects as go
import plotly.express as px
import polars as pl
from RNApysoforms.utils import check_df
from typing import List

def make_transcript_expression_traces(
    long_expression_matrix: pl.DataFrame,
    y: str = "transcript_id",
    x: str = "counts"
) -> list:
    
    boxplot_traces = []

    unique_y

    for sp in species:
        trace = dict(
            y=long_expression_matrix.select(pl.col(y)),
            name=sp,
            boxmean=True,
            marker_color=colors[sp],
            line=dict(width=2)
        ))