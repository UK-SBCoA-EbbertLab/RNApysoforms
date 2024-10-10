import polars as pl
import os

def read_gtf(path: str) -> pl.DataFrame:
    """
    Reads a GTF (Gene Transfer Format) file and returns the data as a Polars DataFrame.

    This function parses a GTF file to extract genomic features, specifically focusing on
    'exon' and 'CDS' (Coding DNA Sequence) feature types. It extracts key attributes from the
    'attributes' column, such as `gene_id`, `gene_name`, `transcript_id`, `transcript_name`,
    `transcript_biotype`, and `exon_number`. The function performs several validation checks
    on the file path and file content, and handles missing values by substituting default values
    where necessary. The resulting DataFrame includes detailed information about each transcript
    feature, suitable for downstream genomic analyses.

    **Expected Columns in GTF File:**
    The GTF file is expected to have no header and the following tab-separated columns:
        1. `seqnames` (chromosome or sequence name)
        2. `source` (annotation source)
        3. `type` (feature type, e.g., 'exon', 'CDS')
        4. `start` (start position of the feature)
        5. `end` (end position of the feature)
        6. `score` (score value, often '.')
        7. `strand` (strand information, '+' or '-')
        8. `phase` (reading frame phase)
        9. `attributes` (semicolon-separated key-value pairs)

    Parameters
    ----------
    path : str
        The file path to the GTF file to be read. The file must have a '.gtf' extension.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame containing extracted gene and transcript features. The DataFrame includes
        the following columns:
        - `gene_id`: Identifier for the gene.
        - `gene_name`: Name of the gene. If missing, filled with `gene_id`.
        - `transcript_id`: Identifier for the transcript.
        - `transcript_name`: Name of the transcript. If missing, filled with `transcript_id`.
        - `transcript_biotype`: Biotype classification of the transcript.
        - `seqnames`: Chromosome or sequence name.
        - `strand`: Strand information ('+' or '-').
        - `type`: Feature type ('exon' or 'CDS').
        - `start`: Start position of the feature.
        - `end`: End position of the feature.
        - `exon_number`: Exon number within the transcript, cast to Int64.

    Raises
    ------
    ValueError
        If the file path does not exist.
        If the path is not a file.
        If the file does not have a '.gtf' extension.
        If required columns are missing or the file cannot be read properly.

    Examples
    --------
    Read a GTF file and display the first few rows:

    >>> from RNApysoforms.io import read_gtf
    >>> df = read_gtf("/path/to/file.gtf")
    >>> print(df.head())

    Notes
    -----
    - The function uses lazy evaluation for reading and processing the file, which is efficient for large GTF files.
    - It filters out feature types other than 'exon' and 'CDS' to focus on relevant transcript features.
    - Regular expressions are used to extract specific attributes from the 'attributes' column.
    - Missing `gene_name` and `transcript_name` values are filled with `gene_id` and `transcript_id`, respectively.
    - The 'exon_number' field is cast to Int64, handling possible nulls without strict type enforcement.
    - The function returns a collected Polars DataFrame after all lazy operations are executed.

    """

    # Validate the file path to ensure it exists and is a file
    if not os.path.exists(path):
        raise ValueError(f"File '{path}' does not exist. Please provide a valid file path.")

    if not os.path.isfile(path):
        raise ValueError(f"'{path}' is not a file. Please provide a valid file path.")

    # Ensure the file has a '.gtf' extension
    if not path.lower().endswith('.gtf'):
        raise ValueError("File must have a '.gtf' extension.")

    # Define the column names and data types for the GTF file
    column_names = [
        "seqnames",    # Chromosome or sequence name
        "source",      # Annotation source
        "type",        # Feature type (e.g., exon, CDS)
        "start",       # Start position of the feature
        "end",         # End position of the feature
        "score",       # Score value (usually '.')
        "strand",      # Strand information ('+' or '-')
        "phase",       # Reading frame phase
        "attributes"   # Additional attributes in key-value pairs
    ]
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

    # Lazily read the GTF file using Polars for efficiency with large files
    lazy_df = pl.scan_csv(
        path,
        separator="\t",
        has_header=False,
        comment_prefix="#",       # Skip comment lines starting with '#'
        new_columns=column_names, # Assign column names since GTF files have no header
        dtypes=dtypes             # Specify data types for each column
    )

    # Filter for features of interest: 'exon' and 'CDS'
    lazy_filtered = lazy_df.filter(pl.col("type").is_in(["exon", "CDS"]))

    # Extract attributes from the 'attributes' column using regular expressions
    lazy_extracted = lazy_filtered.with_columns([
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id"),
        pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),
        pl.col("attributes").str.extract(r'transcript_name "([^"]+)"', 1).alias("transcript_name"),
        pl.col("attributes").str.extract(r'transcript_biotype "([^"]+)"', 1).alias("transcript_biotype"),
        pl.col("attributes").str.extract(r'exon_number "([^"]+)"', 1).alias("exon_number")
    ])

    # Fill missing 'gene_name' and 'transcript_name' with 'gene_id' and 'transcript_id' respectively
    lazy_filled = lazy_extracted.with_columns([
        pl.col("gene_name").fill_null(pl.col("gene_id")),
        pl.col("transcript_name").fill_null(pl.col("transcript_id"))
    ])

    # Select and reorder the relevant columns for the final DataFrame
    result_lazy = lazy_filled.select([
        "gene_id",
        "gene_name",
        "transcript_id",
        "transcript_name",
        "transcript_biotype",
        "seqnames",
        "strand",
        "type",
        "start",
        "end",
        "exon_number"
    ])

    # Cast 'exon_number' to Int64, handling possible nulls without strict type enforcement
    result_lazy = result_lazy.with_columns([
        pl.col("exon_number").cast(pl.Int64, strict=False)
    ])

    # Collect the DataFrame to execute all lazy operations and get the final result
    result_df = result_lazy.collect()

    return result_df  # Return the processed Polars DataFrame
