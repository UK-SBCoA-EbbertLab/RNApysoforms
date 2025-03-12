import polars as pl
from RNApysoforms.utils import check_df
from RNApysoforms.calculate_exon_number import calculate_exon_number

def to_intron(annotation: pl.DataFrame, transcript_id_column: str = "transcript_id") -> pl.DataFrame:
    """
    Converts exon coordinates into corresponding intron coordinates within a genomic annotation dataset.

    This function identifies introns by calculating the genomic intervals between consecutive exons for each transcript.
    It returns a DataFrame with the calculated intron coordinates and retains relevant grouping based on the specified
    `transcript_id_column`, typically 'transcript_id'. If intron entries are already present
     in the input data, the function exits with an error.

    Parameters
    ----------
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations, including exon coordinates with the following required columns:
        - `seqnames`: Chromosome or sequence name.
        - `start`: Start position of the exon.
        - `end`: End position of the exon.
        - `type`: Feature type, expected to include "exon".
        - `exon_number`: Numerical identifier for exons.
        - `transcript_id_column`: Column used to group transcripts, typically "transcript_id".
    transcript_id_column : str, optional
        The column used to group data, typically 'transcript_id'. Default is 'transcript_id'.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame containing both exon and intron coordinates, including other genomic features such as CDS if present.
        The DataFrame includes the following columns:
        - `seqnames`
        - `start`
        - `end`
        - `type` ("exon" or "intron")
        - `exon_number`
        - Additional columns from the input DataFrame.

    Raises
    ------
    TypeError
        If `annotation` is not a Polars DataFrame.
    ValueError
        If the input DataFrame does not contain the required columns (`seqnames`, `start`, `end`, `type`, `exon_number`, and `transcript_id_column`).
        If the input DataFrame already contains introns.

    Examples
    --------
    Convert exons into introns:

    >>> import polars as pl
    >>> from RNApysoforms import to_intron
    >>> df = pl.DataFrame({
    ...    "seqnames": ["chr1", "chr1", "chr1"],
    ...    "start": [100, 200, 300],
    ...    "end": [150, 250, 350],
    ...    "type": ["exon", "exon", "exon"],
    ...    "transcript_id": ["tx1", "tx1", "tx1"],
    ...    "strand": ["+", "+", "+"],
    ...    "exon_number": [1, 2, 3]
    ... })

    >>> df_with_introns = to_intron(df, transcript_id_column="transcript_id")
    >>> print(df_with_introns.head())

    This will return a DataFrame with calculated intron positions between the provided exon coordinates.

    Notes
    -----
    - The function filters out invalid introns where `start` or `end` is null, and introns with length ≤ 1 are discarded.
    - The input DataFrame must contain the required columns listed above.
    - If 'exon_number' is not present in the input DataFrame, it will be automatically calculated.
    - The function can handle input DataFrames with or without existing intron entries. If intron entries are absent, the function generates them.
    - Additional genomic features (e.g., CDS) present in the input DataFrame are retained and merged with intron entries.
    - The function does not adjust intron positions by adding or subtracting 1; intron positions are directly taken from exon boundaries.

    """

    # Check if annotation is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}."
            "\nYou can convert a pandas DataFrame to Polars using: polars_df = pl.from_pandas(pandas_df)"
        )

    # Validate the input DataFrame to ensure required columns are present
    check_df(annotation, ["seqnames", "start", "strand", "end", "type", transcript_id_column])

    if "exon_number" not in annotation.columns:
        annotation = calculate_exon_number(annotation, transcript_id_column)

    ## Make sure annotation has no introns
    if not annotation.filter(pl.col("type") == "intron").is_empty():
        raise ValueError("Your annotation already has introns, please get rid of them before using this function")
    
    ## Make sure annotation has no introns
    if annotation.filter(pl.col("type") == "exon").is_empty():
        raise ValueError("Your annotation must contains exon entries.")
    
    # Separate exons from the annotation data
    exons = annotation.filter(pl.col("type") == "exon").clone()

    # Sort exons by transcript ID and start position
    exons_sorted = exons.sort([transcript_id_column, "start"])

    # Add the end position of the previous exon within each transcript
    exons_with_shift = exons_sorted.with_columns([
        pl.col("end").shift(1).over(transcript_id_column).alias("prev_end"),
        pl.col("exon_number").shift(1).over(transcript_id_column).alias("prev_exon_number"),
    ])

    # Check for overlapping exons
    overlaps = exons_with_shift.filter(
        (pl.col("start") <= pl.col("prev_end")) &
        (pl.col("prev_end").is_not_null())
    ).select([
        transcript_id_column, "exon_number", "start", "end",
        pl.col("prev_exon_number"), pl.col("prev_end")
    ])

    # If overlaps are found, raise an error
    if not overlaps.is_empty():
        raise ValueError(
            f"Your annotation has exons with overlapping coordinates for the same {transcript_id_column} value. "
            "Two exons from the same transcript should not overlap; your annotation may be corrupted. "
            f"Here are the problematic entries:\n{overlaps}"
        )


    # Separate exons and other features (e.g., CDS) from the annotation data
    exons = annotation.filter(pl.col("type") == "exon")
    other_features = annotation.filter(pl.col("type") != "exon")

    # Sort exons by transcript ID and genomic coordinates to ensure correct intron calculation
    sort_columns = [transcript_id_column, 'start', 'end']
    exons_sorted = exons.sort(sort_columns)

    # Calculate intron start and end positions by shifting exon coordinates within each transcript group
    exons_with_introns = exons_sorted.with_columns([
        (pl.col('end').shift(1).over(transcript_id_column) + 1).alias('intron_start'),  # Intron start = end of previous exon + 1 (GTF coordinates)
        (pl.col('start') - 1).alias('intron_end'),                                      # Intron end = start of current exon - 1 (GTF coordinates)
        pl.col("exon_number").shift(1).over(transcript_id_column).alias('intron_number'), ## Get intron number
        pl.lit('intron').alias('type')                                            # Set feature type as 'intron'
    ])

    # Exclude columns that are either renamed or already processed
    exclude_cols = ['start', 'end', 'intron_start', 'intron_end', 'type', 'exon_number']
    columns_to_add = [col for col in exons.columns if col not in exclude_cols]

    # Handle additional columns by taking the first value in each group if transcript_id_column exists
    if transcript_id_column:
        other_cols_expr = [pl.col(col).first().over(transcript_id_column).alias(col) for col in columns_to_add]
    else:
        other_cols_expr = [pl.col(col).first().alias(col) for col in columns_to_add]

    # Select intron columns and include any additional required columns
    introns = exons_with_introns.select([
        pl.col('intron_start').alias('start'),  # Intron start position
        pl.col('intron_end').alias('end'),      # Intron end position
        pl.col("intron_number").alias("exon_number"),                  # Retain exon_number column for reference
        pl.col('type'),                         # Type of feature ('intron')
        *other_cols_expr                        # Include additional columns as necessary
    ])

    # Fix exon number for negative strand introns
    introns = introns.with_columns(
    pl.when(pl.col("strand") == "-")
    .then(pl.col("exon_number") - 1)
    .otherwise(pl.col("exon_number"))
    .alias("exon_number")
    )


    # Remove rows where either 'start' or 'end' is null (invalid introns)
    introns = introns.drop_nulls(subset=['start', 'end'])

    # Filter out introns where the length is 1 or less (invalid introns)
    introns = introns.filter((pl.col('end') - pl.col('start')).abs() > 1)

    # Cast 'start' and 'end' columns to integers for genomic coordinates
    introns = introns.with_columns([
        pl.col('start').cast(pl.Int64),
        pl.col('end').cast(pl.Int64)
    ])

    # Reorder intron columns to match the order of exons for consistency
    introns = introns[exons.columns]

    # Concatenate exons, other features, and introns into a single DataFrame
    combined_annotation = pl.concat([exons, other_features, introns])

    # Sort the combined DataFrame by 'seqnames', transcript_id_column, 'start', 'end', and 'type' for organized output
    combined_annotation = combined_annotation.sort(
        ["seqnames", transcript_id_column, "start", "end", "type"], descending=False
    )

    return combined_annotation  # Return the combined DataFrame with intron entries
