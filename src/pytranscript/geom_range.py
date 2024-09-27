# import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
# import pandas as pd  # Imports pandas for data manipulation and analysis
# import numpy as np  # Imports numpy for numerical operations

# def geom_range(data, x_start='start', x_end='end', y='transcript_name',
#                fill='grey', color='black', opacity=1, line_width=0.25, height=0.5):
#     """
#     Create a list of Plotly range (exon) trace objects.
    
#     Parameters:
#     -----------
#     data : pd.DataFrame
#         DataFrame containing the range (exon) data.
#     x_start : str, optional
#         Column name for the start position of exons. Default is 'start'.
#     x_end : str, optional
#         Column name for the end position of exons. Default is 'end'.
#     y : str, optional
#         Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
#     fill : str or list, optional
#         Fill color or a list of colors. Default is 'grey'.
#     color : str, optional
#         Color of the rectangle borders. Default is 'black'.
#     opacity : float, optional
#         Opacity of the rectangles. Default is 1.
#     line_width : float, optional
#         Width of the rectangle borders. Default is 0.25.
#     height : float, optional
#         Height of the rectangles. Default is 0.5.
    
#     Returns:
#     --------
#     list
#         A list of Plotly trace objects representing the exons.
#     """
#     traces = []  # Initialize an empty list to store Plotly trace objects (rectangles representing exons)
    
#     # Get the unique values from the 'y' column (e.g., transcript names) to assign them distinct y-positions
#     unique_y = data[y].unique()  
    
#     # Create a dictionary to map each unique transcript name to a numerical y-position
#     y_dict = {val: i for i, val in enumerate(unique_y)}
    
#     # If 'fill' is a single string, convert it into a list of the same color for all data points
#     if isinstance(fill, str):
#         fill = [fill] * len(data)  
    
#     # Iterate over each row in the DataFrame to create a Plotly rectangle for each exon
#     for idx, row in data.iterrows():
#         y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript
#         # Define the rectangle trace (an exon) with position and appearance attributes
#         trace = dict(
#             type="rect",  # Specifies that this trace is a rectangle
#             x0=row[x_start],  # Start position of the exon (x-axis, left boundary)
#             x1=row[x_end],  # End position of the exon (x-axis, right boundary)
#             y0=y_pos - height/2,  # Bottom boundary of the rectangle (y-axis)
#             y1=y_pos + height/2,  # Top boundary of the rectangle (y-axis)
#             fillcolor=fill[idx],  # The color used to fill the rectangle (exon)
#             line=dict(color=color, width=line_width),  # Border color and width of the rectangle
#             opacity=opacity  # Transparency of the rectangle
#         )
#         # Append the created trace to the 'traces' list
#         traces.append(trace)
    
#     # Return the list of Plotly trace objects
#     return traces


### TEST CASE THAT WORKS BY FIXING  GENOMIC COORDINATES ###

# Sample DataFrame representing exons with start and end positions, and transcript names
#data = pd.DataFrame({
#    'start': [100, 200, 300],
#    'end': [150, 250, 350],
#    'transcript_name': ['Transcript_1', 'Transcript_2', 'Transcript_1']
#})

# Run the function to generate the traces
#traces = geom_range(data, x_start='start', x_end='end', y='transcript_name')

# Create a Plotly figure
#fig = go.Figure()

# Add each trace (rectangle) to the figure
#for trace in traces:
#    fig.add_shape(trace)  # Adding the rectangles using 'add_shape' function

# Define the range for x-axis to match exon start and end values
#x_min = data['start'].min() - 10  # Padding for better view
#x_max = data['end'].max() + 10  # Padding for better view

# Customize the layout of the plot with correct x and y range
#fig.update_layout(
#    title="Exon Ranges",
#    xaxis_title="Genomic Position",
#    yaxis_title="Transcripts",
#    yaxis=dict(
#        tickmode='array',
#        tickvals=[0, 1],  # Corresponding y positions for transcript names
#        ticktext=['Transcript_1', 'Transcript_2']  # Labels for the y positions
#    ),
#    xaxis=dict(
#        range=[x_min, x_max],  # Set the x-axis range explicitly based on exon start/end
#        tickmode='linear',  # Ensure linear scaling on x-axis
#        dtick=50  # Set tick intervals as per the genomic data scale (adjust if needed)
#    ),
#    height=400,  # Adjust height to fit the number of transcripts
#    showlegend=False
#)

# Display the figure
#fig.show()

import plotly.graph_objects as go  # Imports Plotly for creating interactive plots
import polars as pl  # Imports Polars for data manipulation and analysis

def geom_range(
    data,
    x_start="start",
    x_end="end",
    y="transcript_name",
    fill="grey",
    color="black",
    opacity=1,
    line_width=0.25,
    height=0.5,
):
    """
    Create a list of Plotly range (exon) trace objects.

    Parameters:
    -----------
    data : pl.DataFrame
        DataFrame containing the range (exon) data.
    x_start : str, optional
        Column name for the start position of exons. Default is 'start'.
    x_end : str, optional
        Column name for the end position of exons. Default is 'end'.
    y : str, optional
        Column name for the y-axis values (e.g., transcript names). Default is 'transcript_name'.
    fill : str or list, optional
        Fill color or a list of colors. Default is 'grey'.
    color : str, optional
        Color of the rectangle borders. Default is 'black'.
    opacity : float, optional
        Opacity of the rectangles. Default is 1.
    line_width : float, optional
        Width of the rectangle borders. Default is 0.25.
    height : float, optional
        Height of the rectangles. Default is 0.5.

    Returns:
    --------
    list
        A list of Plotly trace objects representing the exons.
    """
    traces = []  # Initialize an empty list to store Plotly trace objects (rectangles representing exons)

    # Get the unique values from the 'y' column (e.g., transcript names) to assign them distinct y-positions
    unique_y = data[y].unique(maintain_order=True).to_list()

    # Create a dictionary to map each unique transcript name to a numerical y-position
    y_dict = {val: i for i, val in enumerate(unique_y)}

    # If 'fill' is a single string, convert it into a list of the same color for all data points
    if isinstance(fill, str):
        fill = [fill] * len(data)

    # Iterate over each row in the DataFrame to create a Plotly rectangle for each exon
    for idx, row in enumerate(data.iter_rows(named=True)):
        y_pos = y_dict[row[y]]  # Get the corresponding y-position for the current transcript
        # Define the rectangle trace (an exon) with position and appearance attributes
        trace = dict(
            type="rect",  # Specifies that this trace is a rectangle
            x0=row[x_start],  # Start position of the exon (x-axis, left boundary)
            x1=row[x_end],  # End position of the exon (x-axis, right boundary)
            y0=y_pos - height / 2,  # Bottom boundary of the rectangle (y-axis)
            y1=y_pos + height / 2,  # Top boundary of the rectangle (y-axis)
            fillcolor=fill[idx],  # The color used to fill the rectangle (exon)
            line=dict(color=color, width=line_width),  # Border color and width of the rectangle
            opacity=opacity,  # Transparency of the rectangle
        )
        # Append the created trace to the 'traces' list
        traces.append(trace)

    # Return the list of Plotly trace objects
    return traces