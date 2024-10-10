import polars as pl
from typing import List

def check_df(df: pl.DataFrame, required_cols: List[str]):
    """
    Validates that the input Polars DataFrame contains all the required columns.

    This utility function checks whether the provided Polars DataFrame includes all columns specified in the
    `required_cols` list. It is commonly used to ensure that DataFrames meet the necessary schema requirements
    before processing.

    Parameters
    ----------
    df : pl.DataFrame
        The Polars DataFrame to validate.
    required_cols : List[str]
        A list of required column names that the DataFrame must contain.

    Raises
    ------
    TypeError
        If `df` is not a Polars DataFrame.
    ValueError
        If any of the required columns are missing from the DataFrame.

    Examples
    --------
    Validate a DataFrame for required columns:

    >>> import polars as pl
    >>> from RNApysoforms.utils import check_df
    >>> df = pl.DataFrame({
    ...     "col1": [1, 2, 3],
    ...     "col2": ["a", "b", "c"]
    ... })
    >>> check_df(df, required_cols=["col1", "col2"])
    # No exception is raised; all required columns are present.

    >>> check_df(df, required_cols=["col1", "col3"])
    Traceback (most recent call last):
    ValueError: The DataFrame is missing the following required columns: col3

    Notes
    -----
    - The function first checks if `df` is an instance of `pl.DataFrame`.
    - It then compares the columns in `df` against the `required_cols` list.
    - If any required columns are missing, it raises a `ValueError` listing the missing columns.
    - This function is useful for input validation in data processing pipelines.
    - No return value is provided; the function raises exceptions if validation fails.

    """
    
    # Ensure the input is a Polars DataFrame
    if not isinstance(df, pl.DataFrame):
        raise TypeError(
            f"Expected 'df' to be of type pl.DataFrame, got {type(df)}."
            "\nYou can convert a pandas DataFrame to Polars using: polars_df = pl.from_pandas(pandas_df)"
        )

    # Identify any missing columns by comparing against the DataFrame's columns
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    # Raise an error if there are missing columns
    if missing_cols:
        raise ValueError(
            f"The DataFrame is missing the following required columns: {', '.join(missing_cols)}"
        )
