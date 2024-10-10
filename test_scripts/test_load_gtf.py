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
                                       order_by_expression=True, keep_top_expressed_transcripts=5)

#MIR99AHG
#SOD1
#RUNX1

## Shorten gaps
rescaled_annotation = pt.shorten_gaps(annotation=annotation, transcript_id_column="transcript_id")


## Create traces
traces = pt.make_traces(
    expression_matrix = counts,
    order_transcripts_by_expression_matrix=True,
    y='transcript_id', 
    expression_columns=["counts", "relative_abundance"],
    x_start="rescaled_start",
    x_end="rescaled_end",
    annotation_fill_color="blue",
    expression_fill_color="grey",
    exon_height=0.3,
    cds_height=0.5,
    annotation_hue="transcript_biotype",
    expression_hue=None,
    arrow_height=0.2,
    arrow_length=1.6,
    expression_plot_style="boxplot", 
    spanmode="hard")


fig = pt.make_plot(traces = traces,
                    subplot_titles = ["Transcript Structure", "Counts"],
                    showlegend = True,
                    height = 900,
                    width = 1800,
                    boxgap=0.2,
                    boxgroupgap=0,
                    horz_grid_expression_plot=False)





# Show or save the plot
fig.show()