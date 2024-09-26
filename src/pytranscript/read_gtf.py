import polars as pl
import re
import os

def _extract_attributes(attribute_str: str) -> dict:
    """
    Extracts gene_id, gene_name, transcript_id, and transcript_name from the attribute string in GTF format.
    
    The attributes column in a GTF file contains multiple key-value pairs such as:
    gene_id "value"; gene_name "value"; transcript_id "value"; transcript_name "value";
    
    This function uses regular expressions to extract these values, if present.

    Args:
        attribute_str (str): The attributes field from a GTF file row.

    Returns:
        dict: A dictionary containing the extracted 'gene_id', 'gene_name', 'transcript_id', and 'transcript_name'.
              If a specific attribute is not present, its value will be None.
    """
    # Search for each attribute in the format 'key "value"' using regular expressions
    gene_id = re.search(r'gene_id "([^"]+)"', attribute_str)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_str)
    transcript_id = re.search(r'transcript_id "([^"]+)"', attribute_str)
    transcript_name = re.search(r'transcript_name "([^"]+)"', attribute_str)
    
    # Return the attributes as a dictionary, using None if any attribute is not found
    return {
        'gene_id': gene_id.group(1) if gene_id else None,
        'gene_name': gene_name.group(1) if gene_name else None,
        'transcript_id': transcript_id.group(1) if transcript_id else None,
        'transcript_name': transcript_name.group(1) if transcript_name else None
    }

def read_gtf(path: str) -> pl.DataFrame:
    """
    Reads a GTF file, filters for entries of type "exon" or "CDS", extracts necessary attributes,
    and returns a structured dataframe with selected columns.

    The function reads a GTF file, which is a tab-delimited file format used for representing gene structures.
    It filters the entries to only include those with types "exon" and "CDS", then extracts key attributes
    such as gene_id, gene_name, transcript_id, and transcript_name from the attributes column.
    
    Additionally, if gene_name is missing (NA), it will be filled with gene_id. Similarly, if transcript_name
    is missing (NA), it will be filled with transcript_id.

    Args:
        path (str): The file path to the GTF file. The file must have a .gtf extension.

    Returns:
        pl.DataFrame: A polars DataFrame with the columns:
                      ["gene_id", "gene_name", "transcript_id", "transcript_name", "chr", "strand", "start", "end"]
    
    Raises:
        ValueError: If the file path is invalid, the file extension is not .gtf, or the file cannot be read.
    """
    
    # Check if the file path exists
    if not os.path.exists(path):
        raise ValueError(f"File '{path}' does not exist. Please provide a valid file path.")
    
    # Check if the provided path is a file
    if not os.path.isfile(path):
        raise ValueError(f"'{path}' is not a file. Please provide a valid file path.")
    
    # Check if the file has a .gtf extension
    if not path.endswith('.gtf'):
        raise ValueError("File must have a .gtf extension.")
    
    # Try reading the GTF file using polars
    try:
        df = pl.read_csv(
            path, 
            sep="\t",  # GTF files are tab-separated
            has_header=False,  # GTF files don't have headers
            comment_char="#",  # Lines beginning with # are comments and should be ignored
            new_columns=["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]  # Define the columns
        )
    except Exception as e:
        raise ValueError(f"Failed to read the GTF file: {e}")
    
    # Filter the dataframe to keep only rows where the 'type' column is either "exon" or "CDS"
    df_filtered = df.filter(pl.col("type").is_in(["exon", "CDS"]))
    
    # Apply the attribute extraction function to the 'attributes' column for each row
    # This extracts gene_id, gene_name, transcript_id, and transcript_name for each entry
    attributes_extracted = df_filtered["attributes"].apply(lambda attr: _extract_attributes(attr))
    
    # Convert the list of extracted dictionaries into a DataFrame with columns for each attribute
    extracted_df = pl.DataFrame(attributes_extracted.tolist())
    
    # Add the extracted attributes as columns to the filtered dataframe
    final_df = df_filtered.with_columns([
        extracted_df["gene_id"], 
        extracted_df["gene_name"], 
        extracted_df["transcript_id"], 
        extracted_df["transcript_name"]
    ])
    
    # Fill NAs in the 'gene_name' column with values from 'gene_id'
    final_df = final_df.with_columns(
        pl.when(pl.col("gene_name").is_null())
          .then(pl.col("gene_id"))
          .otherwise(pl.col("gene_name"))
          .alias("gene_name")
    )
    
    # Fill NAs in the 'transcript_name' column with values from 'transcript_id'
    final_df = final_df.with_columns(
        pl.when(pl.col("transcript_name").is_null())
          .then(pl.col("transcript_id"))
          .otherwise(pl.col("transcript_name"))
          .alias("transcript_name")
    )
    
    # Select and return the final set of columns in the specified order
    final_df = final_df.select([
        "gene_id", "gene_name", "transcript_id", "transcript_name", "chr", "strand", "start", "end"
    ])
    
    return final_df
