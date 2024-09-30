import polars as pl
from typing import Union, List


def check_df(df: pl.DataFrame, required_cols: List[str]):
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
    if not isinstance(df, pl.DataFrame):
        raise ValueError("Object must be a Polars DataFrame. "
                            "Other object types are currently not supported. "
                            "If you are using a Pandas DataFrame, please convert to Polars before proceding.")


    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"DataFrame must have columns: {', '.join(missing_cols)}")