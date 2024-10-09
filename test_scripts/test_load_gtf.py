## Import libraries
import RNApysoforms as pt
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


## Read gtf
annotation = pt.read_gtf("./raw_data/Homo_sapiens.GRCh38.110.gtf")


counts = pt.read_expression_matrix(expression_matrix_path="./test_data/counts_matrix_chr21_and_Y.tsv", 
                               metadata_path="./test_data/sample_metadata.tsv",
                               cpm_normalization=False, relative_abundance=True)


# Define a mapping from transcript_biotype to colors
biotype_colors = {
    'protein_coding': '#F8766D',
    'retained_intron': '#00BFC4',
    'protein_coding_CDS_not_defined': 'green'
}


## Define gene name to filter
gene_name = "APP"

## Filter gene name in annotation and counts matrix
annotation, counts = pt.gene_filtering(annotation=annotation, expression_matrix=counts, 
                                       target_gene=gene_name, transcript_id_column="transcript_id",
                                       order_by_expression=False)

print(annotation.head())
print(counts.head())


#MIR99AHG
#SOD1
#RUNX1

## Shorten gaps
rescaled_annotation = pt.shorten_gaps(annotation=annotation, transcript_id_column="transcript_id")

print(rescaled_annotation.head())
print(rescaled_annotation.columns)



## Create traces
traces = pt.make_transcript_structure_traces(
    annotation=rescaled_annotation,
    y='transcript_id', 
    x_start="rescaled_start",
    x_end="rescaled_end",
    fill_color="grey",
    exon_height=0.3,
    cds_height=0.5,
    hue="transcript_biotype",
    arrow_height=0.2,
    arrow_length=1.6,
    
)

## Set expression traces
expression_traces = pt.make_transcript_expression_traces(counts, hue="AD status", expression_columns=["counts", "relative_abundance"])

fig = pt.make_plot(transcript_structure_traces = traces, expression_traces=expression_traces,
                   annotation=rescaled_annotation, 
                    subplot_titles = ["Transcript Structure", "Counts", "Relative Abundance"],
                    showlegend = True,
                    height = 900,
                    width = 1800)


# Update layout and show the plot
# fig.update_layout(
#     row=1,
#     col=1,
#     title={'text': f"{gene_name} Transcript Structure", 'x': 0.5, 'y': 0.9, 'xanchor': 'center', 'yanchor': 'top',
#            'font': dict(family='DejaVu Sans', size=14)},
#     xaxis_title="",
#     yaxis_title="",
#     height=400,
#     width=800,
#     showlegend=True,
#     yaxis=dict(
#         tickvals=list(range(counts.n_unique(subset="transcript_id"))),               # Positions of the ticks
#         ticktext=counts.select("transcript_id").unique(maintain_order=True).to_series().to_list(),  # Custom labels for the ticks
#         tickfont=dict(size=10, family='DejaVu Sans', color='black')),
#     xaxis=dict(showticklabels=False),
#     legend=dict(traceorder="reversed"),
#     hoverlabel=dict(
#     font=dict(
#         size=8  # Set the desired hover font size here
#     )
#     )
#     )



# Update layout and show the plot
# fig.update_layout(
#     row=1,
#     col=2,
#     title={'text': f"{gene_name} Counts", 'x': 0.5, 'y': 0.9, 'xanchor': 'center', 'yanchor': 'top',
#            'font': dict(family='DejaVu Sans', size=14)},
#     xaxis_title="",
#     yaxis_title="",
#     height=400,
#     width=800,
#     showlegend=True,
#     yaxis=dict(
#         tickvals=list(range(rescaled_annotation.n_unique(subset="transcript_id"))),               # Positions of the ticks
#         ticktext=rescaled_annotation.select("transcript_id").unique(maintain_order=True).to_series().to_list(),  # Custom labels for the ticks
#         tickfont=dict(size=10, family='DejaVu Sans', color='black')),
#     xaxis=dict(showticklabels=False),
#     legend=dict(traceorder="reversed"),
#     hoverlabel=dict(
#     font=dict(
#         size=8  # Set the desired hover font size here
#     )
#     )
#     )


# Show or save the plot
fig.show()

#fig.write_html(("transcript_structure.html"))