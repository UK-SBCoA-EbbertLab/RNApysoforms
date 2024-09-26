import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def geom_intron(data, x_start='start', x_end='end', y='transcript_name', 
                strand='strand', color='black', line_width=0.5, opacity=1,
                arrow_color='black', arrow_size=0.3, arrow_min_intron_length=100):
    """
    Create a list of Plotly intron trace objects with strand arrows.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the intron data.
    x_start : str, optional
        Column name for the start position of introns. Default is 'start'.
    x_end : str, optional
        Column name for the end position of introns. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    strand : str, optional
        Column name for the strand information. Default is 'strand'.
    color : str, optional
        Color of the intron lines. Default is 'black'.
    line_width : float, optional
        Width of the intron lines. Default is 0.5.
    opacity : float, optional
        Opacity of the intron lines. Default is 1.
    arrow_color : str, optional
        Color of the strand arrows. Default is 'black'.
    arrow_size : float, optional
        Size of the strand arrows. Default is 0.3.
    arrow_min_intron_length : int, optional
        Minimum intron length to draw an arrow. Default is 0.
    
    Returns:
    --------
    list
        A list of Plotly trace objects representing the introns with optional strand arrows.
    """
    traces = []
    unique_y = data[y].unique()
    y_dict = {val: i for i, val in enumerate(unique_y)}
    
    for idx, row in data.iterrows():
        y_pos = y_dict[row[y]]
        trace = dict(
            type="line",
            x0=row[x_start],
            x1=row[x_end],
            y0=y_pos,
            y1=y_pos,
            line=dict(color=color, width=line_width),
            opacity=opacity
        )
        traces.append(trace)

        ## Define arrow direction
        arrow_direction = ">" if row['strand'] == "+" else "<"

        if (row[x_end] - row[x_start]) > arrow_min_intron_length:

            if arrow_min_intron_length < 1:
                arrow_min_intron_length=100

            num_arrows = int((row[x_end] - row[x_start]) / (arrow_min_intron_length))
            for i in range(num_arrows):
                arrow_x = (row[x_start] + i * ((row[x_end] - row[x_start]) / num_arrows))
                if ((abs(arrow_x - row[x_end]) > 100) and (abs(arrow_x - row[x_start]) > 100)):
                    arrow_trace = go.Scatter(
                        x=[arrow_x],
                        y=[y_pos],
                        mode='markers',
                        marker=dict(symbol="arrow-right" if row["strand"] == "+" else "arrow-left",  # Use 'symbol' as the key
                        size=(arrow_size * 10),
                        color="black"),
                        showlegend=False,
                        hoverinfo="none"
                    )
                    traces.append(arrow_trace)
    
    return traces