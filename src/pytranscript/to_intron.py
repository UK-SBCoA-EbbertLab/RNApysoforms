import polars as pl
from pytranscript.utils import check_df  # Utility function for DataFrame validation

def to_intron(annotation: pl.DataFrame, group_var: str = "transcript_id") -> pl.DataFrame:

    """
    Converts exon coordinates into corresponding intron coordinates.

    This function identifies the introns between exons in a genomic annotation dataset by shifting exon coordinates. 
    It returns a DataFrame with the calculated intron positions, maintaining any relevant groupings like transcript IDs.

    Parameters:
        exons (pl.DataFrame): A DataFrame containing exon coordinates with 'start' and 'end' columns.
        group_var (str, optional): Column used to group transcripts, default is 'transcript_id'.

    Returns:
        pl.DataFrame: A DataFrame containing intron coordinates derived from the exon data.
    """

    # Check required columns in the input DataFrame
    if group_var is None:
        check_df(annotation, ["seqnames", "start", "end", "type"])
    else:
        check_df(annotation, ["seqnames", "start", "end", "type", group_var])

    ## Separate CDS and exon
    exons = annotation.filter(pl.col("type") == "exon")
    cds = annotation.filter(pl.col("type") == "CDS")


    # Sort exons by group_var, start, and end coordinates
    sort_cols =  [group_var, 'start', 'end']
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
    annotation = pl.concat([exons, cds, introns])

    ## Sort
    annotation = annotation.sort(["seqnames", group_var, "start", "end", "type"], descending=False)

    return annotation  # Return the DataFrame annotation containing valid intron coordinates

