import polars as pl
from typing import Union, List

def to_intron(exons: pl.DataFrame, group_var: Union[str, List[str]] = None) -> pl.DataFrame:
    """
    Convert exon coordinates to intron coordinates.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon coordinates. Must have 'start' and 'end' columns.
    group_var : str or List[str], optional
        Column name(s) for grouping transcripts. Default is None.

    Returns:
    --------
    pl.DataFrame
        DataFrame containing the intron coordinates.
    """
    
    ## Define required columns
    required_cols = ["start", "end"]

    # Check if all required columns are present in the exons DataFrame
    assert all(col in exons.columns for col in required_cols), \
        f"exons DataFrame must have columns: {', '.join(required_cols)}"

    if group_var is not None:
        if isinstance(group_var, str):
            group_var = [group_var]
        if not all(col in exons.columns for col in group_var):
            raise ValueError("group_var columns must exist in exons DataFrame")
    else:
        group_var = []

    # Sort exons by start and end coordinates
    sort_cols = group_var + ['start', 'end']
    exons_sorted = exons.sort(sort_cols)

    # Calculate intron start, end, and intron number
    exons_with_introns = exons_sorted.with_columns([
        pl.col('end').shift(1).over(group_var).alias('intron_start'),
        pl.col('start').alias('intron_end'),
        pl.lit('intron').alias('type'),
        pl.col("exon_number").alias("exon_number")
    ])

    # Exclude columns that have been renamed or already included
    exclude_cols = ['start', 'end', 'intron_start', 'intron_end', 'type', 'exon_number']
    columns_to_add = [col for col in exons.columns if col not in exclude_cols]
    

    # For other columns, get the first value per group
    if group_var:
        other_cols_expr = [pl.col(col).first().over(group_var).alias(col) for col in columns_to_add]
    else:
        other_cols_expr = [pl.col(col).first().alias(col) for col in columns_to_add]


    # Select necessary columns with unique aliases to avoid duplication
    introns = exons_with_introns.select([
        pl.col('intron_start').alias('start'),
        pl.col('intron_end').alias('end'),
        pl.col("exon_number"),
        pl.col('type'),
        *[expr.alias(f'{expr.meta.root_names()[0]}') for i, expr in enumerate(other_cols_expr)]  # Ensure unique aliases for other expressions
    ])

    # Remove NAs and introns with width of 1
    introns = introns.drop_nulls(subset=['start', 'end'])
    introns = introns.filter((pl.col('end') - pl.col('start')).abs() != 1)

    # Cast from float to integer
    introns = introns.with_columns([
        pl.col('start').cast(pl.Int64),
        pl.col('end').cast(pl.Int64)
    ])

    return introns
