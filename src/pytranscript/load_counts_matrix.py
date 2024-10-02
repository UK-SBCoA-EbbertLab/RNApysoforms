import polars as pl
from typing import List, Union, Optional
import warnings

def load_counts_matrix(
    counts_path: str,
    metadata_path: Optional[str] = None,
    cpm_normalization: bool = False,
    counts_feature_id_columns: Union[List[str], str] = ["gene_id", "transcript_id"],
    metadata_sample_id_column: str = "sample_id"
) -> pl.DataFrame:
   
    """
    Load and process a counts matrix, optionally merging metadata and performing CPM normalization.

    This function reads a counts matrix file and, optionally, a metadata file, merging the two on a sample identifier. 
    The counts can also be normalized to Counts Per Million (CPM) for further analysis. The resulting DataFrame is returned
    in long format, with counts and optional CPM values, merged with metadata if provided.

    Parameters
    ----------
    counts_path : str
        Path to the counts matrix file.
    metadata_path : str, optional
        Path to the metadata file. If provided, the metadata will be merged with the counts data.
    cpm_normalization : bool, optional
        Whether to perform Counts Per Million (CPM) normalization on the counts data. Default is False.
    counts_feature_id_columns : list of str or str, optional
        Column name(s) in the counts DataFrame that identify features (e.g., genes or transcripts).
        Default is ["gene_id", "transcript_id"].
    metadata_sample_id_column : str, optional
        Column name in the metadata DataFrame that identifies samples. Default is "sample_id".

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame in long format with the counts data (and CPM values if normalization is performed),
        optionally merged with metadata.

    Raises
    ------
    ValueError
        If required columns are missing or if data types are incorrect.
    Warning
        If there is partial overlap between samples in counts data and metadata.

    Examples
    --------
    Load a counts matrix and perform CPM normalization, merging with metadata:

    >>> df = load_counts_matrix("counts.csv", metadata_path="metadata.csv", cpm_normalization=True)
    >>> print(df.head())

    Load a counts matrix without normalization:

    >>> df = load_counts_matrix("counts.csv")
    >>> print(df.head())

    Notes
    -----
    - The function assumes that the counts data and metadata are stored in CSV format and will read them accordingly.
    - If CPM normalization is performed, the counts will be scaled to reflect Counts Per Million for each feature across samples.
    - The function will raise warnings if there is partial sample overlap between counts data and metadata.
    - The resulting DataFrame will be in long format, with counts or CPM data for each sample-feature combination.
    """

    # Load counts_path as a polars dataframe.
    try:
        counts_df = pl.read_csv(counts_path)
    except Exception as e:
        raise ValueError(f"Failed to read counts data from '{counts_path}': {e}")

    # Ensure counts_feature_id_columns is a list
    if isinstance(counts_feature_id_columns, str):
        counts_feature_id_columns = [counts_feature_id_columns]

    # Check if counts_feature_id_columns are present in counts_dataframe
    missing_columns = [col for col in counts_feature_id_columns if col not in counts_df.columns]
    if missing_columns:
        raise ValueError(f"The following feature ID columns are missing in the counts dataframe: {missing_columns}")

    # Check if all columns that are not counts_feature_id_columns are numerical
    counts_columns = [col for col in counts_df.columns if col not in counts_feature_id_columns]
    non_numeric_columns = [col for col in counts_columns if not pl.datatypes.is_numeric_dtype(counts_df[col].dtype)]
    if non_numeric_columns:
        raise ValueError(f"The following columns are expected to be numerical but are not: {non_numeric_columns}")

    if cpm_normalization:
        # Perform CPM normalization for each transcript
        counts_df = counts_df.with_columns([
            (
                pl.col(col) / pl.col(col).sum() * 1e6
            ).alias(col + "_CPM")
            for col in counts_columns
        ])

        # Transform dataframe into long format
        # Melt counts columns
        counts_long = counts_df.melt(
            id_vars=counts_feature_id_columns,
            value_vars=counts_columns,
            variable_name=metadata_sample_id_column,
            value_name="counts"
        )

        # Melt CPM columns
        CPM_columns = [col + "_CPM" for col in counts_columns]
        CPM_long = counts_df.melt(
            id_vars=counts_feature_id_columns,
            value_vars=CPM_columns,
            variable_name=metadata_sample_id_column,
            value_name="CPM"
        ).with_column(
            pl.col(metadata_sample_id_column).str.replace(r"_CPM$", "")
        )

        # Join counts_long and CPM_long
        long_counts_df = counts_long.join(
            CPM_long,
            on=counts_feature_id_columns + [metadata_sample_id_column],
            how="left"
        )

    else:
        # Transform dataframe into long format
        long_counts_df = counts_df.melt(
            id_vars=counts_feature_id_columns,
            value_vars=counts_columns,
            variable_name=metadata_sample_id_column,
            value_name="counts"
        )

    if metadata_path is not None:
        # Load metadata_path as a polars dataframe
        try:
            metadata_df = pl.read_csv(metadata_path)
        except Exception as e:
            raise ValueError(f"Failed to read metadata from '{metadata_path}': {e}")

        # Check if metadata_sample_id_column is present in metadata_dataframe
        if metadata_sample_id_column not in metadata_df.columns:
            raise ValueError(f"The metadata_sample_id_column '{metadata_sample_id_column}' is not present in the metadata dataframe.")

        # Check if any values in the metadata_sample_id_column are present in the long_counts_dataframe metadata_sample_id_column
        counts_sample_ids = long_counts_df[metadata_sample_id_column].unique().to_list()
        metadata_sample_ids = metadata_df[metadata_sample_id_column].unique().to_list()

        overlapping_sample_ids = set(counts_sample_ids).intersection(set(metadata_sample_ids))

        if not overlapping_sample_ids:
            raise ValueError("No overlapping sample IDs found between counts data and metadata.")

        # Check how many values in the metadata_sample_id_column are not present in the counts data
        metadata_sample_ids_not_in_counts = set(metadata_sample_ids) - set(counts_sample_ids)
        counts_sample_ids_not_in_metadata = set(counts_sample_ids) - set(metadata_sample_ids)

        warning_message = ""
        if metadata_sample_ids_not_in_counts:
            warning_message += f"The following sample IDs are present in metadata but not in counts data: {list(metadata_sample_ids_not_in_counts)}. "
        if counts_sample_ids_not_in_metadata:
            warning_message += f"The following sample IDs are present in counts data but not in metadata: {list(counts_sample_ids_not_in_metadata)}."
        if warning_message:
            warnings.warn(warning_message)

        # Merge metadata_dataframe to long_counts_dataframe
        long_counts_df = long_counts_df.join(
            metadata_df,
            on=metadata_sample_id_column,
            how="left"
        )

    return long_counts_df
