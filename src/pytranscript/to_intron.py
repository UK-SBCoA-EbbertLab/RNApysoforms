import polars as pl
from pytranscript.utils import check_df  # Utility function for DataFrame validation

def to_intron(exons: pl.DataFrame, group_var: str = "transcript_id") -> pl.DataFrame:

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
        check_df(exons, ["start", "end"])
    else:
        check_df(exons, ["start", "end", group_var])

    # Define the required columns for exons
    required_cols = ["start", "end"]
    
    # Ensure required columns are present in the exons DataFrame
    assert all(col in exons.columns for col in required_cols), \
        f"exons DataFrame must have columns: {', '.join(required_cols)}"

    # If group_var is specified, ensure it's a list, and validate its columns
    if group_var is not None:
        if isinstance(group_var, str):
            group_var = [group_var]  # Convert string to list for consistent processing
        if not all(col in exons.columns for col in group_var):
            raise ValueError("group_var columns must exist in exons DataFrame")
    else:
        group_var = []  # Set group_var as an empty list if not provided

    # Sort exons by group_var, start, and end coordinates
    sort_cols = group_var + ['start', 'end']
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
        pl.col("exon_number"),  # Exon number (if applicable)
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

    return introns  # Return the DataFrame containing valid intron coordinates

