## Import libraries
import RNApysoforms as pt
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


## Read gtf
annotation = pt.read_gtf("./raw_data/Homo_sapiens.GRCh38.110.gtf")


counts = pt.read_counts_matrix(counts_path="./test_data/counts_matrix_chr21_and_Y.tsv", 
                               metadata_path="./test_data/sample_metadata.tsv",
                               cpm_normalization=True)

# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': '#F8766D',
    'retained_intron': '#00BFC4',
    'protein_coding_CDS_not_defined': 'green'
}


## Add biotype colors
annotation = annotation.with_columns(pl.col('transcript_biotype').replace_strict(biotype_colors, default="gray").alias('fillcolor'))


## Define gene name to filter
gene_name = "APP"

## Filter gene name in annotation and counts matrix
annotation, counts = pt.gene_filtering(annotation=annotation, counts_matrix=counts, 
                                       gene_name_to_filter=gene_name, group_var="transcript_id")


#MIR99AHG
#SOD1
#RUNX1

## Shorten gaps
rescaled_annotation = pt.shorten_gaps(annotation=annotation, group_var="transcript_id")

## Create traces
traces = pt.make_traces(
    data=rescaled_annotation,
    y='transcript_id', 
    x_start="rescaled_start",
    x_end="rescaled_end",
    fill_color="grey",
    exon_height=0.3,
    cds_height=0.5,
    hue="transcript_biotype",
    is_hoverable=True,
    arrow_height=0.5,
    arrow_length=1,
    #color_palette=px.colors.qualitative.Plotly,
    #color_map=biotype_colors
)


# Create the plot
fig = go.Figure()

# Create the plot and add all traces in one step for efficiency
fig = go.Figure(data=traces)


# Call the new function to set the genomic axis range
fig = pt.set_axis(fig, data=rescaled_annotation, x_start="rescaled_start", x_end="rescaled_end")


# Update layout and show the plot
fig.update_layout(
    title={'text': f"{gene_name} Transcript Structure", 'x': 0.5, 'y': 0.9, 'xanchor': 'center', 'yanchor': 'top',
           'font': dict(family='DejaVu Sans', size=14)},
    xaxis_title="",
    yaxis_title="",
    height=400,
    width=800,
    showlegend=True,
    yaxis=dict(
        tickvals=list(range(rescaled_annotation.n_unique(subset="transcript_id"))),               # Positions of the ticks
        ticktext=rescaled_annotation.select("transcript_id").unique(maintain_order=True).to_series().to_list(),  # Custom labels for the ticks
        tickfont=dict(size=10, family='DejaVu Sans', color='black')),
    xaxis=dict(showticklabels=False),
    legend=dict(traceorder="reversed"),
    hoverlabel=dict(
    font=dict(
        size=8  # Set the desired hover font size here
    )
    )
    )

# Show or save the plot
fig.show()
#fig.write_html(("transcript_structure.html"))