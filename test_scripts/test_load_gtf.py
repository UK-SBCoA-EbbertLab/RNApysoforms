## Import libraries
import pytranscript as pt
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## Read gtf
annotations = pt.read_gtf("./test_data/Homo_sapiens.GRCh38.112.chr21-22.gtf")


# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': '#F8766D',
    'retained_intron': '#00BFC4',
    'protein_coding_CDS_not_defined': 'green'
}

## Add biotype colors
annotations = annotations.with_columns(pl.col('transcript_biotype').replace_strict(biotype_colors, default="gray").alias('fillcolor'))

## Define gene
gene_name = "RUNX1"

## Filter gene of interest
annotations = annotations.filter(pl.col("gene_name") == gene_name)

#MIR99AHG

## Shorten gaps
rescaled_annotations = pt.shorten_gaps(annotation=annotations, group_var="transcript_id")

## Create traces
traces = pt.make_traces(
    data=rescaled_annotations,
    x_start='start',
    x_end='end',
    y='transcript_id', 
    fill_column="fillcolor",
    fill_color="grey",
    exon_height=0.3,
    cds_height=0.5,
    strand="strand"
)

# Create the plot
fig = go.Figure()

#Batch add exon, CDS, and intron shapes
fig.update_layout(
    shapes=[trace for trace in traces if isinstance(trace, dict)]
)

# Call the new function to set the genomic axis range
fig = pt.set_axis(fig, data=rescaled_annotations)


# Update layout and show the plot
fig.update_layout(
    title={'text': f"{gene_name} Transcript Structure", 'x': 0.5, 'y': 0.9, 'xanchor': 'center', 'yanchor': 'top',
           'font': dict(family='DejaVu Sans', size=14)},
    xaxis_title="",
    yaxis_title="",
    height=400,
    width=800,
    showlegend=False,
    yaxis=dict(
        tickvals=list(range(rescaled_annotations.n_unique(subset="transcript_id"))),               # Positions of the ticks
        ticktext=rescaled_annotations.select("transcript_id").unique(maintain_order=True).to_series().to_list(),  # Custom labels for the ticks
        tickfont=dict(size=10, family='DejaVu Sans', color='black')),
        xaxis=dict(showticklabels=False)
)

# Show or save the plot
fig.show()
fig.write_html(("transcript_structure.html"))