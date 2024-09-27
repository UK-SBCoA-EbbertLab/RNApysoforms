import polars as pl
from typing import Union, List

def check_coord_object(df: pl.DataFrame, check_seqnames: bool = False, check_strand: bool = False):
    """
    Check if the DataFrame has the required columns for coordinate data.

    Parameters:
    -----------
    df : pl.DataFrame
        The DataFrame to check.
    check_seqnames : bool, optional
        Whether to check for the 'seqnames' column. Default is False.
    check_strand : bool, optional
        Whether to check for the 'strand' column. Default is False.

    Raises:
    -------
    ValueError
        If the DataFrame doesn't meet the required structure.
    """
    if not isinstance(df, pl.DataFrame):
        raise ValueError("Object must be a Polars DataFrame. "
                         "Other object types are currently not supported.")

    required_cols = ['start', 'end']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        raise ValueError(f"DataFrame must have the columns: {', '.join(missing_cols)}")
    
    if check_seqnames and 'seqnames' not in df.columns:
        raise ValueError("DataFrame must have the column 'seqnames'")
    
    if check_strand and 'strand' not in df.columns:
        raise ValueError("DataFrame must have the column 'strand'")

def check_group_var(df: pl.DataFrame, group_var: Union[str, List[str]]):
    """
    Check if the group_var exists in the DataFrame.

    Parameters:
    -----------
    df : pl.DataFrame
        The DataFrame to check.
    group_var : str or List[str]
        The grouping variable(s) to check for.

    Raises:
    -------
    ValueError
        If any of the group_var columns are missing from the DataFrame.
    """
    if group_var:
        if isinstance(group_var, str):
            group_var = [group_var]
        missing_vars = [var for var in group_var if var not in df.columns]
        if missing_vars:
            raise ValueError(f"group_var columns missing from DataFrame: {', '.join(missing_vars)}")

def check_required_columns(df: pl.DataFrame, required_cols: List[str], context: str = "DataFrame"):
    """
    Check if the DataFrame has all the required columns.

    Parameters:
    -----------
    df : pl.DataFrame
        The DataFrame to check.
    required_cols : List[str]
        List of required column names.
    context : str, optional
        Context for the error message. Default is "DataFrame".

    Raises:
    -------
    ValueError
        If any of the required columns are missing.
    """
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"{context} must have columns: {', '.join(missing_cols)}")

def get_unique_values(df: pl.DataFrame, column: str) -> List:
    """
    Get unique values from a DataFrame column.

    Parameters:
    -----------
    df : pl.DataFrame
        The DataFrame to check.
    column : str
        The column name to get unique values from.

    Returns:
    --------
    List
        List of unique values in the column.
    """
    return df[column].unique().to_list()

def calculate_width(df: pl.DataFrame) -> pl.Series:
    """
    Calculate the width of genomic features.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing 'start' and 'end' columns.

    Returns:
    --------
    pl.Series
        Series containing the calculated widths.
    """
    return df['end'] - df['start'] + 1

def sort_genomic_coordinates(df: pl.DataFrame, group_var: Union[str, List[str]] = None) -> pl.DataFrame:
    """
    Sort DataFrame by genomic coordinates.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame to sort.
    group_var : str or List[str], optional
        Column name(s) for grouping before sorting.

    Returns:
    --------
    pl.DataFrame
        Sorted DataFrame.
    """
    if group_var is None:
        sort_cols = ['start', 'end']
    else:
        if isinstance(group_var, str):
            group_var = [group_var]
        sort_cols = group_var + ['start', 'end']
    return df.sort(by=sort_cols)



# Sample DataFrame
df = pl.DataFrame({
    'seqnames': ['chr1', 'chr1', 'chr2'],
    'start': [100, 150, 200],
    'end': [500, 550, 600],
    'strand': ['+', '-', '+'],
    'gene': ['gene1', 'gene2', 'gene3']
})

# Test check_coord_object
check_coord_object(df, check_seqnames=True, check_strand=True)

# Test check_group_var
check_group_var(df, 'gene')

# Test check_required_columns
check_required_columns(df, ['start', 'end', 'gene'])

# Test get_unique_values
unique_genes = get_unique_values(df, 'gene')
print(unique_genes)  # Output: ['gene1', 'gene2', 'gene3']

# Test calculate_width
widths = calculate_width(df)
print(widths)  # Output: pl.Series with calculated widths

# Test sort_genomic_coordinates
sorted_df = sort_genomic_coordinates(df, group_var='seqnames')
print(sorted_df)