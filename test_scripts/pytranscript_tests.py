import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import functions from our plotly_ggtranscript package
import sys
import os
sys.path.append("C:/Users/local_bag222/Desktop/dash_apps/AD_RNAseq_dash_app")
import pytranscript

def get_valid_input(df, column_name):
    """
    Continuously prompts the user for input until a valid value (one that exists in the specified column of the DataFrame) is provided.
    
    Parameters:
    df (pd.DataFrame): The DataFrame in which to check for the input value.
    column_name (str): The column name of the DataFrame where the value will be searched.
    
    Returns:
    str: The valid input provided by the user that exists in the DataFrame column.
    """
    
    # Infinite loop to keep asking for input until a valid one is provided
    while True:
        # Request input from the user
        user_input = input(f"Enter a value for {column_name}: ")
        
        # Check if the input exists in the specified column of the DataFrame
        if user_input in df[column_name].values:
            print(f"'{user_input}' found in the column '{column_name}'.")
            return user_input  # Return the valid input
        
        # If the input does not exist in the column, ask the user to try again
        else:
            print(f"'{user_input}' not found in column '{column_name}'. Please try again.")




## Import annotations
annotation = pd.read_csv("./raw_data/test_data/test_annotations.tsv", sep="\t")

print(annotation.head(10))

## Get gene name
gene_name = get_valid_input(annotation, "gene_name")
annotation = annotation.loc[annotation["gene_name"] == gene_name].copy()

print(annotation.head())


# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': '#F8766D',
    'processed_transcript': '#00BFC4',
    'NA': 'gray'
}



## Add biotype colors
annotation['fillcolor'] = annotation['transcript_biotype'].map(biotype_colors)

# Extract exons
exons = annotation[annotation['type'] == 'exon']

## Extract CDS
cds = annotation[annotation['type'] == 'CDS']

## Calculate CDS and exon differences
cds_diff = calculate_cds_exon_difference(cds, exons)


# Use .loc to avoid SettingWithCopyWarning
exons = exons.copy()

## Rescale SOD exons
rescaled = shorten_gaps(exons=exons, introns=to_intron(exons, "transcript_name"), group_var = "transcript_name")

## Define rescaled exons and introns
rescaled_exons = rescaled.loc[rescaled["type"] == "exon"].copy()
rescaled_introns = rescaled.loc[rescaled["type"] == "intron"].copy()

## Correct CDS coordinates for new exon coordinates
rescaled_cds = rescale_cds(cds_diff, rescaled_exons)

# Create the plot
fig = go.Figure()

rescaled_cds.sort_values(by=["transcript_name"], inplace=True)
rescaled_exons.sort_values(by=["transcript_name"], inplace=True)

# Add exons using geom_range, passing the fillcolor directly
exon_traces = geom_range(
    data=rescaled_exons,
    x_start='start',
    x_end='end',
    y='transcript_name',
    fill=rescaled_exons['fillcolor'], 
    height=0.3
)

## Add CDS traces
cds_traces = geom_range(
    data=rescaled_cds,
    x_start='start',
    x_end='end',
    y='transcript_name',
    fill=rescaled_cds['fillcolor'],
    height= 0.5
)


# Create introns and add them using geom_intron
#sod1_introns = to_intron(sod1_exons, group_var="transcript_name")
intron_traces = geom_intron(
    data=rescaled_introns,
    x_start='start',
    x_end='end',
    y='transcript_name',
    strand='strand',
    arrow_min_intron_length=400,
    arrow_size=1
)

# Add exons, CDS, and introns as before
for trace in exon_traces:
    fig.add_shape(trace)

for trace in cds_traces:
    fig.add_shape(trace)

for trace in intron_traces:
    if isinstance(trace, dict):
        fig.add_shape(trace)
    else:
        fig.add_trace(trace)

# Call the new function to set the genomic axis range
fig = set_axis(fig, rescaled_exons, rescaled_introns)

# Update layout and show the plot
fig.update_layout(
    title={'text': f"{gene_name} Transcript Structure", 'x': 0.5, 'y': 0.8, 'xanchor': 'center', 'yanchor': 'top',
           'font': dict(family='DejaVu Sans', size=14)},
    xaxis_title="",
    yaxis_title="",
    height=400,
    width=800,
    showlegend=False,
    yaxis=dict(
        tickvals=list(range(rescaled_exons["transcript_name"].nunique())),               # Positions of the ticks
        ticktext=rescaled_exons["transcript_name"].unique().tolist(),  # Custom labels for the ticks
        tickfont=dict(size=10, family='DejaVu Sans', color='black')),
        xaxis=dict(showticklabels=False)
)

# Show or save the plot
fig.show()
fig.write_html((gene_name + "_transcript_structure.html"))
