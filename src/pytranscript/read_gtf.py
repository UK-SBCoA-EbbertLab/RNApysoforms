import polars as pl  # Polars is used for efficient data manipulation and analysis
import os  # Used to validate file paths

def read_gtf(path: str) -> pl.DataFrame:
   
    """
    Reads a GTF (Gene Transfer Format) file, extracts relevant transcript information, and returns it as a Polars DataFrame.

    The function parses GTF entries, focusing on exon and CDS features, and extracts specific gene and transcript attributes.
    It ensures the GTF file is valid and structured, providing an efficient and structured way to work with transcript data.

    Parameters:
        path (str): Path to the GTF file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing gene and transcript details, such as gene_id, transcript_id, exon_number, and more.

    Raises:
        ValueError: If the file path is invalid or the file is not a GTF.
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
        separator="\t",  # The GTF file is tab-separated
        has_header=False,  # GTF files do not have a header row
        comment_prefix="#",  # Ignore comment lines that start with '#'
        new_columns=["seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        # Map the GTF columns to new names for easier reference
    )
    
    # Filter the LazyFrame to keep rows where the 'type' column is either 'exon' or 'CDS'
    lazy_filtered = lazy_df.filter(pl.col("type").is_in(["exon", "CDS"]))
    
    # Extract the required attributes from the 'attributes' column using regular expressions
    lazy_extracted = lazy_filtered.with_columns([
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id"),  # Extract gene_id
        pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),  # Extract gene_name
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),  # Extract transcript_id
        pl.col("attributes").str.extract(r'transcript_name "([^"]+)"', 1).alias("transcript_name"),  # Extract transcript_name
        pl.col("attributes").str.extract(r'transcript_biotype "([^"]+)"', 1).alias("transcript_biotype"),  # Extract transcript_biotype
        pl.col("attributes").str.extract(r'exon_number "([^"]+)"', 1).alias("exon_number")  # Extract exon number
    ])
    
    # Fill missing gene_name values with gene_id, and missing transcript_name with transcript_id
    lazy_final = lazy_extracted.with_columns([
        pl.when(pl.col("gene_name").is_null())  # If gene_name is missing
          .then(pl.col("gene_id"))  # Use gene_id as a fallback
          .otherwise(pl.col("gene_name"))  # Otherwise keep the original gene_name
          .alias("gene_name"),
        
        pl.when(pl.col("transcript_name").is_null())  # If transcript_name is missing
          .then(pl.col("transcript_id"))  # Use transcript_id as a fallback
          .otherwise(pl.col("transcript_name"))  # Otherwise keep the original transcript_name
          .alias("transcript_name")
    ])
    
    # Select the final columns to include in the output DataFrame
    result_lazy = lazy_final.select([
        "gene_id", "gene_name", "transcript_id", "transcript_name", 
        "transcript_biotype", "seqnames", "strand", "type", "start", "end", "exon_number"
    ])

    # Ensure the start, end, and exon_number columns are stored as integers
    result_lazy = result_lazy.with_columns([
        pl.col("start").cast(pl.Int64),  # Convert start coordinates to integer
        pl.col("end").cast(pl.Int64),  # Convert end coordinates to integer
        pl.col("exon_number").cast(pl.Float64)  # Convert exon_number to integer
    ])
    
    # Trigger the lazy computation by collecting the data into an eager DataFrame
    return result_lazy.collect()

