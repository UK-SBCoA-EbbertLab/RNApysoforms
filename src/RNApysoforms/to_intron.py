import polars as pl
from RNApysoforms.utils import check_df

def to_intron(annotation: pl.DataFrame, group_var: str = "transcript_id") -> pl.DataFrame:
   
    """
    Converts exon coordinates into corresponding intron coordinates within a genomic annotation dataset.

    This function identifies introns by calculating the genomic intervals between consecutive exons for each transcript.
    It returns a DataFrame with the calculated intron coordinates and retains relevant grouping based on the specified 
    `group_var`, typically 'transcript_id'.

    Parameters
    ----------
    annotation : pl.DataFrame
        A DataFrame containing genomic annotations, including exon coordinates with group_var, 'start', 'end', and 'exon_number' columns.
    group_var : str, optional
        The column used to group data, typically 'transcript_id'. Default is 'transcript_id'.

    Returns
    -------
    pl.DataFrame
        A DataFrame containing intron coordinates derived from the exons data, as well as other genomic features like CDS.

    Raises
    ------
    ValueError
        If the input DataFrame does not contain the required columns (`seqnames`, `start`, `end`, `type`, `exon_number`, 
        and `group_var`).

    Examples
    --------
    Convert exons into introns:

    >>> import polars as pl
    >>> from RNApysoforms.utils import check_df
    >>> df = pl.DataFrame({
    ...     "seqnames": ["chr1", "chr1", "chr1", "chr1"],
    ...     "start": [100, 200, 400, 600],
    ...     "end": [150, 250, 450, 650],
    ...     "type": ["exon", "exon", "exon", "exon"],
    ...     "transcript_id": ["tx1", "tx1", "tx1", "tx1"],
    ...     "exon_number": [1, 2, 3, 4]
    ... })
    >>> df_with_introns = to_intron(df, group_var="transcript_id")
    >>> print(df_with_introns.head())
    
    This will return a DataFrame with calculated intron positions between the provided exon coordinates.

    Notes
    -----
    - The function filters out invalid introns where `start` or `end` is null, and introns with length ≤ 1 are discarded.
    - The input DataFrame must contain the `seqnames`, `start`, `end`, `type`, `exon_number`, and group_var
    - The resulting DataFrame retains all necessary columns from the input exons and supports optional columns.
    - The input DataFrame can contain just exons or other "type" values as well such as CDS.
    """

    # Check required columns in the input DataFrame
    check_df(annotation, ["seqnames", "start", "end", "type", "exon_number", group_var])

    ## Separate CDS and exon
    exons = annotation.filter(pl.col("type") == "exon")
    other = annotation.filter(pl.col("type") != "exon")

    # Sort exons by group_var, start, and end coordinates
    sort_cols = [group_var, 'start', 'end']
    exons_sorted = exons.sort(sort_cols)

    # Calculate intron start and end positions, shift 'end' to get 'intron_start'
    exons_with_introns = exons_sorted.with_columns([
        pl.col('end').shift(1).over(group_var).alias('intron_start'),  # Intron start = previous exon end
        pl.col('start').alias('intron_end'),  # Intron end = next exon start
        pl.lit('intron').alias('type'),  # Set type as 'intron'
        pl.col("exon_number").alias("exon_number")  # Retain exon_number column
    ])

    # Exclude certain columns that are either renamed or already processed
    exclude_cols = ['start', 'end', 'intron_start', 'intron_end', 'type', 'exon_number']
    columns_to_add = [col for col in exons.columns if col not in exclude_cols]

    # Handle additional columns by taking the first value in each group if group_var exists
    if group_var:
        other_cols_expr = [pl.col(col).first().over(group_var).alias(col) for col in columns_to_add]
    else:
        other_cols_expr = [pl.col(col).first().alias(col) for col in columns_to_add]

    # Select intron columns (start, end, exon_number, type) and any other required columns
    introns = exons_with_introns.select([
        pl.col('intron_start').alias('start'),  # Intron start position
        pl.col('intron_end').alias('end'),  # Intron end position
        (pl.col("exon_number") + 0.5),  # Exon number (if applicable)
        pl.col('type'),  # Type of feature (intron)
        *other_cols_expr  # Include additional columns as necessary
    ])

    # Remove rows where either 'start' or 'end' is null (invalid introns)
    introns = introns.drop_nulls(subset=['start', 'end'])

    # Filter out introns where the length is 1 or less (invalid introns)
    introns = introns.filter((pl.col('end') - pl.col('start')).abs() > 1)

    # Cast 'start' and 'end' columns to integers for genomic coordinates
    introns = introns.with_columns([
        pl.col('start').cast(pl.Int64),
        pl.col('end').cast(pl.Int64)
    ])

    ## Set intron column order to match exons and cds
    introns = introns[exons.columns]

    ## Update annotations to include introns and exons
    annotation = pl.concat([exons, other, introns])

    ## Sort to make output more neat
    annotation = annotation.sort(["seqnames", group_var, "start", "end", "type"], descending=False)

    return annotation  # Return the DataFrame annotation containing valid intron coordinates