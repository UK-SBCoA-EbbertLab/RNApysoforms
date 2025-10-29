import polars as pl
import os

def process_ensembl_gtf(gtf_df: pl.DataFrame) -> pl.DataFrame:
    """
    Processes an already-loaded GTF DataFrame to extract and format genomic features.

    This function takes a Polars DataFrame containing GTF data and processes it by:
    - Filtering for 'exon' and 'CDS' feature types
    - Extracting key attributes (gene_id, gene_name, transcript_id, etc.) from the attributes column
    - Filling missing gene_name and transcript_name with gene_id and transcript_id
    - Validating that required attributes are present (consistent with 2024 ENSEMBL GTF format)
    - Casting exon_number to Int64

    This function is useful when you have already loaded your GTF data into a Polars DataFrame and want to process it
    without reading from a file. If you need to read from a file, use :func:`read_ensembl_gtf` instead.

    **Expected DataFrame Format:**
    The input DataFrame should have the following columns (standard GTF format):
        - `seqnames` (chromosome or sequence name)
        - `source` (annotation source)
        - `type` (feature type, e.g., 'exon', 'CDS')
        - `start` (start position of the feature)
        - `end` (end position of the feature)
        - `score` (score value, often '.')
        - `strand` (strand information, '+' or '-')
        - `phase` (reading frame phase)
        - `attributes` (semicolon-separated key-value pairs)

    **Required Attributes in the attributes column:**
    The function expects ENSEMBL-formatted attributes with the following keys:
        - `gene_id` (required)
        - `transcript_id` (required)
        - `gene_name` (optional, will be filled with gene_id if missing)
        - `transcript_name` (optional, will be filled with transcript_id if missing)
        - `transcript_biotype` (required)
        - `exon_number` (optional for exon/CDS features)

    Parameters
    ----------
    gtf_df : pl.DataFrame
        A Polars DataFrame containing GTF data in standard format with the expected columns.

    Returns
    -------
    pl.DataFrame
        A processed Polars DataFrame containing extracted gene and transcript features. The DataFrame includes
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
        If the GTF DataFrame is not consistent with the 2024 ENSEMBL GTF format (missing required attributes).

    Examples
    --------
    Process an already-loaded GTF DataFrame:

    >>> import polars as pl
    >>> from RNApysoforms import process_ensembl_gtf
    >>> 
    >>> # Assume gtf_df is already loaded with standard GTF columns
    >>> # gtf_df = pl.read_csv("file.gtf", separator="\t", has_header=False, ...)
    >>> 
    >>> processed_df = process_ensembl_gtf(gtf_df)
    >>> print(processed_df.head())

    Notes
    -----
    - The function filters out feature types other than 'exon' and 'CDS'.
    - Regular expressions are used to extract specific attributes from the 'attributes' column.
    - Missing `gene_name` and `transcript_name` values are filled with `gene_id` and `transcript_id`, respectively.
    - The 'exon_number' field is cast to Int64, handling possible nulls without strict type enforcement.
    - If required attributes (gene_id, transcript_id, transcript_biotype) are missing, a ValueError is raised.

    See Also
    --------
    read_ensembl_gtf : Read and process GTF data from a file.
    """

    # Filter for features of interest: 'exon' and 'CDS'
    filtered_df = gtf_df.filter(pl.col("type").is_in(["exon", "CDS"]))

    # Extract attributes from the 'attributes' column using regular expressions
    extracted_df = filtered_df.with_columns([
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id"),
        pl.col("attributes").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),
        pl.col("attributes").str.extract(r'transcript_name "([^"]+)"', 1).alias("transcript_name"),
        pl.col("attributes").str.extract(r'transcript_biotype "([^"]+)"', 1).alias("transcript_biotype"),
        pl.col("attributes").str.extract(r'exon_number "([^"]+)"', 1).alias("exon_number")
    ])

    # Fill missing 'gene_name' and 'transcript_name' with 'gene_id' and 'transcript_id' respectively
    filled_df = extracted_df.with_columns([
        pl.col("gene_name").fill_null(pl.col("gene_id")),
        pl.col("transcript_name").fill_null(pl.col("transcript_id"))
    ])

    # Select and reorder the relevant columns for the final DataFrame
    result_df = filled_df.select([
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

    # Check for any null values in the DataFrame
    if result_df.null_count().select(pl.all().sum()).row(0)[0] > 0:
        raise ValueError(
            "This GTF file is not consistent with the 2024 ENSEMBL GTF format. \n"
            "See this vignette with an example on how to handle other GTF formats: \n"
            "https://rna-pysoforms.readthedocs.io/en/latest/examples/10.dealing_with_different_gtf_files.html"
        )

    # Cast 'exon_number' to Int64, handling possible nulls without strict type enforcement
    result_df = result_df.with_columns([
        pl.col("exon_number").cast(pl.Int64, strict=False)
    ])

    return result_df

def read_ensembl_gtf(path: str) -> pl.DataFrame:
    """
    Reads a GTF (Gene Transfer Format) file and returns the data as a Polars DataFrame.

    This function parses an ENSEMBL GTF file to extract genomic features, specifically focusing on
    'exon' and 'CDS' (Coding DNA Sequence) feature types. It extracts key attributes from the
    'attributes' column, such as `gene_id`, `gene_name`, `transcript_id`, `transcript_name`,
    `transcript_biotype`, and `exon_number`. The function performs several validation checks
    on the file path and file content, and handles missing values by substituting default values
    where necessary. The resulting DataFrame includes detailed information about each transcript
    feature, suitable for downstream genomic analyses.

    **Expected Columns in GTF File:**
    The GTF file is expected to be an ENSEMBL GTF file, have no header and the following tab-separated columns:
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
        The file path to the ENSEMBL GTF file to be read. The file must have a '.gtf' extension.

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

    >>> from RNApysoforms import read_ensembl_gtf
    >>> df = read_ensembl_gtf("/path/to/file.gtf")
    >>> print(df.head())

    Notes
    -----
    - The function uses lazy evaluation for reading and processing the file, which is efficient for large GTF files.
    - It filters out feature types other than 'exon' and 'CDS' to focus on relevant transcript features.
    - Regular expressions are used to extract specific attributes from the 'attributes' column.
    - Missing `gene_name` and `transcript_name` values are filled with `gene_id` and `transcript_id`, respectively.
    - The 'exon_number' field is cast to Int64, handling possible nulls without strict type enforcement.
    - The function returns a collected Polars DataFrame after all lazy operations are executed.
    - An example ENSEMBL GTF file only containing data for human chromosomes 21 and Y can be found here: https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/blob/main/tests/test_data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf

    See Also
    --------
    process_ensembl_gtf : Process an already-loaded GTF DataFrame without reading from a file.
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
        schema_overrides=dtypes             # Specify data types for each column
    )

    # Collect the lazy DataFrame
    gtf_df = lazy_df.collect()

    # Process the GTF DataFrame using the process_ensembl_gtf function
    return process_ensembl_gtf(gtf_df)
