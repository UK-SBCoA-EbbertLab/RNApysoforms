import polars as pl

## Open GTF, read chromosome 21 and save it ##

dtypes = {
    "seqnames": pl.Utf8,     
    "source": pl.Utf8,     
    "type": pl.Utf8,         
    "start": pl.Int64,      
    "end": pl.Int64,         
    "score": pl.Utf8,    
    "strand": pl.Utf8,      
    "phase": pl.Utf8,        
    "attributes": pl.Utf8    
}

gtf = pl.read_csv("./test_data/Homo_sapiens.GRCh38.110.gtf", separator="\t", has_header=False,  
    comment_prefix="#",  new_columns=["seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
    schema_overrides=dtypes)

gtf = gtf.filter(pl.col("seqnames").is_in(["21", "Y"]))

gtf.write_csv("./test_data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf", include_header=False, separator="\t")


## Extract transcript ids in chromosome 21 and Y from GTF ##
transcript_ids = gtf.with_columns([pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id")])["transcript_id"].unique()


## Open counts matrix, change column names, remove some samples, filter to only keep transcript ids in chr 21 nad Y
counts = pl.read_csv("./test_data/counts_transcript.txt", separator="\t")

counts.columns = ['transcript_id', 'gene_id', 'sample_1', 'a', 'b', 'c', 'd', 'e', 
                  'sample_4', 'f', 'sample_7', 'h', 'sample_2', 'sample_6', 'sample_3',
                    'm', 'sample_5', 'o', 'sample_8','p', 'q']

counts = counts.drop(['a', 'b', 'c', 'd', 'e', "f", "h", "m", "o", "p", "q"])

count = counts[["transcript_id", "gene_id",	"sample_1", "sample_2", "sample_3", "sample_4",	"sample_5", "sample_6", "sample_7" ,"sample_8"]]

counts = counts.filter(pl.col("transcript_id").is_in(transcript_ids))

counts.write_csv("./test_data/counts_matrix_chr21_and_Y.tsv", separator="\t")


