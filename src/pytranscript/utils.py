import pandas as pd
from typing import Union, List

def check_coord_object(df: pd.DataFrame, check_seqnames: bool = False, check_strand: bool = False):
    """
    Check if the DataFrame has the required columns for coordinate data.

    Parameters:
    -----------
    df : pd.DataFrame
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
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Object must be a pandas DataFrame. "
                         "Other object types are currently not supported.")

    required_cols = ['start', 'end']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        raise ValueError(f"DataFrame must have the columns: {', '.join(missing_cols)}")
    
    if check_seqnames and 'seqnames' not in df.columns:
        raise ValueError("DataFrame must have the column 'seqnames'")
    
    if check_strand and 'strand' not in df.columns:
        raise ValueError("DataFrame must have the column 'strand'")

def check_group_var(df: pd.DataFrame, group_var: Union[str, List[str]]):
    """
    Check if the group_var exists in the DataFrame.

    Parameters:
    -----------
    df : pd.DataFrame
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

# ... [Keep the previously defined utility functions] ...

def check_required_columns(df: pd.DataFrame, required_cols: List[str], context: str = "DataFrame"):
    """
    Check if the DataFrame has all the required columns.

    Parameters:
    -----------
    df : pd.DataFrame
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

def get_unique_values(df: pd.DataFrame, column: str) -> List:
    """
    Get unique values from a DataFrame column.

    Parameters:
    -----------
    df : pd.DataFrame
        The DataFrame to check.
    column : str
        The column name to get unique values from.

    Returns:
    --------
    List
        List of unique values in the column.
    """
    return df[column].unique().tolist()

def calculate_width(df: pd.DataFrame):
    """
    Calculate the width of genomic features.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing 'start' and 'end' columns.

    Returns:
    --------
    pd.Series
        Series containing the calculated widths.
    """
    return df['end'] - df['start'] + 1

def sort_genomic_coordinates(df: pd.DataFrame, group_var: Union[str, List[str]] = None):
    """
    Sort DataFrame by genomic coordinates.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame to sort.
    group_var : str or List[str], optional
        Column name(s) for grouping before sorting.

    Returns:
    --------
    pd.DataFrame
        Sorted DataFrame.
    """
    sort_cols = (group_var if group_var else []) + ['start', 'end']
    return df.sort_values(sort_cols)
