# test_check_df.py

import pytest
import polars as pl
from RNApysoforms.utils import check_df

def test_check_df_valid_input():
    """
    Test check_df with a valid DataFrame containing all required columns.
    """
    df = pl.DataFrame({
        "col1": [1, 2, 3],
        "col2": ["a", "b", "c"]
    })
    required_cols = ["col1", "col2"]
    # Should not raise any exception
    check_df(df, required_cols)

def test_check_df_missing_columns():
    """
    Test check_df when required columns are missing.
    """
    df = pl.DataFrame({
        "col1": [1, 2, 3]
    })
    required_cols = ["col1", "col2"]
    with pytest.raises(ValueError) as exc_info:
        check_df(df, required_cols)
    assert "The DataFrame is missing the following required columns: col2" in str(exc_info.value)

def test_check_df_invalid_input_type():
    """
    Test check_df when input is not a Polars DataFrame.
    """
    df = {
        "col1": [1, 2, 3],
        "col2": ["a", "b", "c"]
    }
    required_cols = ["col1", "col2"]
    with pytest.raises(TypeError) as exc_info:
        check_df(df, required_cols)
    assert "Expected 'df' to be of type pl.DataFrame" in str(exc_info.value)
