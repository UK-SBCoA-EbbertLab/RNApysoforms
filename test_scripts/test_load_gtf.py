## Import libraries
from pytranscript import read_gtf
import polars as pl

## Read gtf
annotations = read_gtf("./test_data/Homo_sapiens.GRCh38.112.chr21-22.gtf")


## Filter gene of interest
annotation = annotations.filter(pl.col("gene_name") == "APP")

