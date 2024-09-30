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

## Filter gene of interest
annotations = annotations.filter(pl.col("gene_name") == "RUNX1")

#MIR99AHG

## Shorten gaps
rescaled_annotations = pt.shorten_gaps(annotation=annotations, group_var="transcript_id")


## Vizualize data
rescaled_annotations.write_excel("~/Desktop/test.xlsx")

## Separate
rescaled_cds = rescaled_annotations.filter(pl.col("type") == "CDS").clone()
rescaled_exons = rescaled_annotations.filter(pl.col("type") == "exon").clone()
rescaled_introns = rescaled_annotations.filter(pl.col("type") == "intron").clone()


# Create the plot
fig = go.Figure()

#rescaled_cds.sort_values(by=["transcript_id"], inplace=True)
#rescaled_exons.sort_values(by=["transcript_id"], inplace=True)
#rescaled_introns.sort_values(by=["transcript_id"], inplace=True)

print("1")
# Add exons using geom_range, passing the fillcolor directly
exon_traces = pt.geom_range(
    data=rescaled_exons,
    x_start='start',
    x_end='end',
    y='transcript_id', 
    fill=rescaled_exons["fillcolor"],
    height=0.3
)
print("2")

## Add CDS traces
cds_traces = pt.geom_range(
    data=rescaled_cds,
    x_start='start',
    x_end='end',
    y='transcript_id',
    fill=rescaled_cds["fillcolor"],
    height= 0.5
)

print("3")
# Create introns and add them using geom_intron
#sod1_introns = to_intron(sod1_exons, group_var="transcript_name")
intron_traces = pt.geom_intron(
    data= rescaled_introns,
    x_start='start',
    x_end='end',
    y='transcript_id',
    strand='strand',
    arrow_min_intron_length=400,
    arrow_size=0.5
)

print("4")

# Batch add exon, CDS, and intron shapes
fig.update_layout(
    shapes=exon_traces + cds_traces + [trace for trace in intron_traces if isinstance(trace, dict)]
)

# For traces (like Scatter traces), add them all in a single step
fig.add_traces([trace for trace in intron_traces if not isinstance(trace, dict)])

print("5")
# Call the new function to set the genomic axis range
fig = pt.set_axis(fig, rescaled_exons, rescaled_introns)
gene_name = "APP"

print("6")
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
        tickvals=list(range(rescaled_exons.n_unique(subset="transcript_id"))),               # Positions of the ticks
        ticktext=rescaled_exons.select("transcript_id").unique().to_series().to_list(),  # Custom labels for the ticks
        tickfont=dict(size=10, family='DejaVu Sans', color='black')),
        xaxis=dict(showticklabels=False)
)

print("7")
# Show or save the plot
fig.show()
fig.write_html(("transcript_structure.html"))