import polars as pl
import os

def read_gtf(path: str) -> pl.DataFrame:
    """
    Reads a GTF file lazily, filters for entries of type 'exon' or 'CDS',
    extracts necessary attributes, and returns a structured eager DataFrame.
    
    Parameters:
    path (str): The file path to the GTF file.

    Returns:
    pl.DataFrame: A polars DataFrame containing selected columns: 
                  'gene_id', 'gene_name', 'transcript_id', 'transcript_name',
                  'chr', 'strand', 'type', 'start', and 'end'.
    
    Raises:
    ValueError: If the provided file path is invalid or if the file is not a GTF file.
    """

    # Validate the file path existence
    if not os.path.exists(path):
        raise ValueError(f"File '{path}' does not exist. Please provide a valid file path.")
    
    # Ensure that the path is a file (not a directory)
    if not os.path.isfile(path):
        raise ValueError(f"'{path}' is not a file. Please provide a valid file path.")
    
    # Ensure that the file has a .gtf extension
    if not path.endswith('.gtf'):
        raise ValueError("File must have a .gtf extension.")
    
    # Lazily scan the GTF file, treating it as a CSV-like file with tab separators
    lazy_df = pl.scan_csv(
        path, 
        separator="\t", 
        has_header=False, 
        comment_prefix="#",  # Ignore comment lines that start with '#'
        new_columns=["seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )
    
    # Filter the LazyFrame to keep rows where the 'type' column is either 'exon' or 'CDS'
    lazy_filtered = lazy_df.filter(pl.col("type").is_in(["exon", "CDS"]))
    
    # Extract the required attributes from the 'attributes' column using regular expressions
    lazy_extracted = lazy_filtered.with_columns([
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id"),  # Extract gene_id
        pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),  # Extract gene_name
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),  # Extract transcript_id
        pl.col("attributes").str.extract(r'transcript_name "([^"]+)"', 1).alias("transcript_name"),  # Extract transcript_name
        pl.col("attributes").str.extract(r'transcript_biotype "([^"]+)"', 1).alias("transcript_biotype")  # Extract transcript_biotype
    ])
    
    # Fill missing gene_name values with gene_id, and missing transcript_name with transcript_id
    lazy_final = lazy_extracted.with_columns([
        pl.when(pl.col("gene_name").is_null())  # If gene_name is null
          .then(pl.col("gene_id"))  # Use gene_id as fallback
          .otherwise(pl.col("gene_name"))  # Otherwise keep gene_name
          .alias("gene_name"),
        
        pl.when(pl.col("transcript_name").is_null())  # If transcript_name is null
          .then(pl.col("transcript_id"))  # Use transcript_id as fallback
          .otherwise(pl.col("transcript_name"))  # Otherwise keep transcript_name
          .alias("transcript_name")
    ])
    
    # Select the final columns to include in the output DataFrame
    result_lazy = lazy_final.select([
        "gene_id", "gene_name", "transcript_id", "transcript_name", 
        "transcript_biotype", "seqnames", "strand", "type", "start", "end"
    ])

    ## Make coordinates integers
    result_lazy = result_lazy.with_columns([
    pl.col("start").cast(pl.Int64),
    pl.col("end").cast(pl.Int64)])
    
    # Trigger the lazy computation by collecting the data into an eager DataFrame
    return result_lazy.collect()
