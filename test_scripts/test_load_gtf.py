## Import libraries
import pytranscript as pt
import polars as pl
import pandas as pd

## Read gtf
annotations = pt.read_gtf("./test_data/Homo_sapiens.GRCh38.112.chr21-22.gtf")


## Filter gene of interest
annotation = annotations.filter(pl.col("gene_name") == "APP")

print(annotation.columns)

## Shorten gaps
rescaled_annotations = pt.shorten_gaps(annotation=annotation, group_var="transcript_id")


print(rescaled_annotations.head())