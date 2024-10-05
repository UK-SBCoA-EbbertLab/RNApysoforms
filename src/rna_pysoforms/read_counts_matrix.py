import polars as pl
from typing import Optional
import warnings
import os

def load_counts_matrix(
    counts_path: str,
    metadata_path: Optional[str] = None,
    cpm_normalization: bool = False,
    gene_id_column_name: Optional[str] = "gene_id",
    transcript_id_column_name: str = "transcript_id",
    metadata_sample_id_column: str = "sample_id"
) -> pl.DataFrame:
    """
    Load and process a counts matrix, optionally merging metadata, performing CPM normalization,
    and calculating relative transcript abundance.

    This function reads a counts matrix file and, optionally, a metadata file, merging the two on a sample identifier.
    The counts can also be normalized to Counts Per Million (CPM) and relative transcript abundance can be calculated
    for further analysis. The resulting DataFrame is returned in long format, with counts and optional CPM values and
    relative abundance, merged with metadata if provided.

    Parameters
    ----------
    counts_path : str
        Path to the counts matrix file. Supported file formats include .csv, .tsv, .txt, .parquet, and .xlsx.
    metadata_path : str, optional
        Path to the metadata file. If provided, the metadata will be merged with the counts data. Supported file formats
        are the same as for `counts_path`. Default is None.
    cpm_normalization : bool, optional
        Whether to perform Counts Per Million (CPM) normalization on the counts data. Default is False.
    gene_id_column_name : str, optional
        The name of the column in the counts DataFrame that contains gene identifiers. This column will remain fixed during data transformation.
        If provided, relative transcript abundance will be calculated. Default is "gene_id". If set to None, the gene identifier will not be used.
    transcript_id_column_name : str
        The name of the column in the counts DataFrame that contains transcript identifiers. This parameter is required and cannot be None.
        This column will remain fixed during data transformation. Default is "transcript_id".
    metadata_sample_id_column : str, optional
        Column name in the metadata DataFrame that identifies samples. This column is used to merge the metadata and counts data.
        Default is "sample_id".

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame in long format with the counts data (and CPM values and relative abundance if calculated),
        optionally merged with metadata.

    Raises
    ------
    ValueError
        - If 'transcript_id_column_name' is None.
        - If required feature ID columns are missing in the counts file.
        - If the file format is unsupported or the file cannot be read.
        - If required columns in the metadata are missing or sample IDs do not overlap between counts and metadata.
    Warning
        - If there is partial overlap between samples in counts data and metadata.

    Examples
    --------
    Load a counts matrix, perform CPM normalization, calculate relative transcript abundance, and merge with metadata:

    >>> df = load_counts_matrix("counts.csv", metadata_path="metadata.csv", cpm_normalization=True)
    >>> print(df.head())

    Load a counts matrix without normalization but calculate relative transcript abundance:

    >>> df = load_counts_matrix("counts.csv")
    >>> print(df.head())

    Load a counts matrix without calculating relative abundance:

    >>> df = load_counts_matrix("counts.csv", gene_id_column_name=None)
    >>> print(df.head())

    Notes
    -----
    - The `transcript_id_column_name` parameter is required and cannot be None.
    - The function supports multiple file formats (.csv, .tsv, .txt, .parquet, .xlsx) for both counts and metadata files.
    - If CPM normalization is performed, the counts will be scaled to reflect Counts Per Million for each sample.
    - If `gene_id_column_name` is provided, relative transcript abundance is calculated as (transcript_counts / total_gene_counts) * 100.
      If the total gene counts are zero, the relative abundance is set to zero to avoid division by zero errors.
    - Warnings are raised if there is partial sample overlap between counts data and metadata.
    - The resulting DataFrame is returned in long format, with counts, CPM values, and relative abundance for each sample-feature combination.
    """

    # Check if transcript_id_column_name is None and raise an error if so
    if transcript_id_column_name is None:
        raise ValueError("The 'transcript_id_column_name' is required and cannot be None.")
    
    # Load counts_path using the helper function
    counts_df = _get_open_file(counts_path)

    # Build counts_feature_id_columns from gene_id_column_name and transcript_id_column_name
    counts_feature_id_columns = [transcript_id_column_name]
    
    if gene_id_column_name is not None:
        counts_feature_id_columns.append(gene_id_column_name)

    # Check if counts_feature_id_columns are present in counts_df
    missing_columns = [col for col in counts_feature_id_columns if col not in counts_df.columns]
    if missing_columns:
        raise ValueError(f"The following feature ID columns are missing in the counts dataframe: {missing_columns}")

    # Determine counts columns (exclude feature ID columns)
    counts_columns = [col for col in counts_df.columns if col not in counts_feature_id_columns]

    # Check that counts_columns are numeric
    non_numeric_columns = [col for col in counts_columns if not counts_df[col].dtype.is_numeric()]
    if non_numeric_columns:
        raise ValueError(f"The following columns are expected to be numerical but are not: {non_numeric_columns}")

    # Calculate relative transcript abundance if gene_id_column_name is provided
    if gene_id_column_name is not None:
        counts_df = counts_df.with_columns([
            (
                pl.when(pl.col(col).sum().over(gene_id_column_name) == 0)
                .then(0)
                .otherwise((pl.col(col) / pl.col(col).sum().over(gene_id_column_name)) * 100)
                .alias(col + "_relative_abundance")
            )
            for col in counts_columns
        ])

    if cpm_normalization:
        # Perform CPM normalization for each sample
        counts_df = counts_df.with_columns([
            (
                pl.col(col) / pl.col(col).sum() * 1e6
            ).alias(col + "_CPM")
            for col in counts_columns
        ])

    # Transform counts_df into long format for counts
    counts_long = counts_df.melt(
        id_vars=counts_feature_id_columns,
        value_vars=counts_columns,
        variable_name=metadata_sample_id_column,
        value_name="counts"
    )

    # Initialize long_counts_df with counts_long
    long_counts_df = counts_long

    # If CPM normalization was performed, melt CPM columns and join
    if cpm_normalization:
        CPM_columns = [col + "_CPM" for col in counts_columns]

        CPM_long = counts_df.melt(
            id_vars=counts_feature_id_columns,
            value_vars=CPM_columns,
            variable_name=metadata_sample_id_column,
            value_name="CPM"
        ).with_columns(
            pl.col(metadata_sample_id_column).str.replace(r"_CPM$", "")
        )

        # Join counts_long and CPM_long
        long_counts_df = long_counts_df.join(
            CPM_long,
            on=counts_feature_id_columns + [metadata_sample_id_column],
            how="left"
        )

    # If relative abundance was calculated, melt and join
    if gene_id_column_name is not None:
        relative_abundance_columns = [col + "_relative_abundance" for col in counts_columns]

        relative_abundance_long = counts_df.melt(
            id_vars=counts_feature_id_columns,
            value_vars=relative_abundance_columns,
            variable_name=metadata_sample_id_column,
            value_name="relative_abundance"
        ).with_columns(
            pl.col(metadata_sample_id_column).str.replace(r"_relative_abundance$", "")
        )

        # Join with long_counts_df
        long_counts_df = long_counts_df.join(
            relative_abundance_long,
            on=counts_feature_id_columns + [metadata_sample_id_column],
            how="left"
        )

    if metadata_path is not None:
        # Load metadata_path using the helper function
        metadata_df = _get_open_file(metadata_path)

        # Check if metadata_sample_id_column is present in metadata_df
        if metadata_sample_id_column not in metadata_df.columns:
            raise ValueError(f"The metadata_sample_id_column '{metadata_sample_id_column}' is not present in the metadata dataframe.")

        # Check overlap of sample IDs between counts data and metadata
        counts_sample_ids = long_counts_df[metadata_sample_id_column].unique().to_list()
        metadata_sample_ids = metadata_df[metadata_sample_id_column].unique().to_list()

        overlapping_sample_ids = set(counts_sample_ids).intersection(set(metadata_sample_ids))

        if not overlapping_sample_ids:
            raise ValueError("No overlapping sample IDs found between counts data and metadata.")

        # Warn about sample ID mismatches
        metadata_sample_ids_not_in_counts = set(metadata_sample_ids) - set(counts_sample_ids)
        counts_sample_ids_not_in_metadata = set(counts_sample_ids) - set(metadata_sample_ids)

        warning_message = ""
        if metadata_sample_ids_not_in_counts:
            warning_message += f"The following sample IDs are present in metadata but not in counts data: {list(metadata_sample_ids_not_in_counts)}. "
        if counts_sample_ids_not_in_metadata:
            warning_message += f"The following sample IDs are present in counts data but not in metadata: {list(counts_sample_ids_not_in_metadata)}."
        if warning_message:
            warnings.warn(warning_message)

        # Merge metadata_df to long_counts_df
        long_counts_df = long_counts_df.join(
            metadata_df,
            on=metadata_sample_id_column,
            how="left"
        )

    return long_counts_df

import polars as pl
import os

def _get_open_file(file_path: str) -> pl.DataFrame:
    """
    Open a file based on its extension and load it into a Polars DataFrame.

    This function supports multiple file formats such as .csv, .tsv, .txt, .parquet, and .xlsx. 
    It automatically determines the correct method to open the file based on its extension.

    Parameters
    ----------
    file_path : str
        The path to the file to be opened.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame containing the contents of the file.

    Raises
    ------
    ValueError
        If the file extension is unsupported or the file cannot be read due to an error.

    Examples
    --------
    Open a CSV file:

    >>> df = _get_open_file("data.csv")
    
    Open a TSV file:

    >>> df = _get_open_file("data.tsv")

    Notes
    -----
    - This function is used internally by `load_counts_matrix` to load counts and metadata files.
    - It handles different file extensions and raises a clear error if the file format is unsupported.
    """
    _, file_extension = os.path.splitext(file_path)
    
    try:
        if file_extension in [".tsv", ".txt"]:
            return pl.read_csv(file_path, separator="\t")
        elif file_extension == ".csv":
            return pl.read_csv(file_path)
        elif file_extension == ".parquet":
            return pl.read_parquet(file_path)
        elif file_extension == ".xlsx":
            return pl.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported file extension '{file_extension}'. Supported extensions are .tsv, .txt, .csv, .parquet, .xlsx")
    except Exception as e:
        raise ValueError(f"Failed to read the file '{file_path}': {e}")
