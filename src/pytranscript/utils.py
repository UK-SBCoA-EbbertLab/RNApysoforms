import polars as pl
from typing import List

def check_df(df: pl.DataFrame, required_cols: List[str]):
    """
    Check if the Polars DataFrame contains all the required columns.

    Parameters:
    -----------
    df : pl.DataFrame
        The Polars DataFrame to check.
    required_cols : List[str]
        A list of required column names that the DataFrame must contain.

    Raises:
    -------
    ValueError
        If the input is not a Polars DataFrame, or if any of the required columns are missing.
    """
    
    # Ensure the input is a Polars DataFrame
    if not isinstance(df, pl.DataFrame):
        raise ValueError(
            "Input must be a Polars DataFrame. "
            "If you're using a Pandas DataFrame, please convert it to Polars."
        )

    # Identify any missing columns by comparing against the DataFrame's columns
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    # Raise an error if there are missing columns
    if missing_cols:
        raise ValueError(
            f"The DataFrame is missing the following required columns: {', '.join(missing_cols)}"
        )
