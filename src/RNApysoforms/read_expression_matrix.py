import polars as pl
from typing import Optional, List
import warnings
import os

def read_expression_matrix(
    expression_matrix_path: str,
    metadata_path: Optional[str] = None,
    expression_measure_name: str = "counts",
    cpm_normalization: bool = False,
    relative_abundance: bool = False,
    gene_id_column_name: Optional[str] = "gene_id",
    transcript_id_column_name: str = "transcript_id",
    metadata_sample_id_column: str = "sample_id"
) -> pl.DataFrame:
    """
    Loads and processes an expression matrix, optionally merging with metadata, performing CPM normalization, and calculating relative transcript abundance.

    This function reads an expression matrix file and, optionally, a metadata file, merging the two on a specified sample identifier column.
    It supports performing Counts Per Million (CPM) normalization and calculating relative transcript abundance based on gene counts.
    The resulting DataFrame is returned in long format, including the expression measures, optional CPM values, relative abundances, and metadata if provided.

    **Expression Matrix Format Requirements:**
    - Must be in wide format with samples as columns and transcripts as rows
    - First column(s) must contain the feature identifiers (transcript_id and optionally gene_id)
    - Sample columns must contain numeric expression values (e.g., read counts)
    - Sample column names must exactly match the values in the metadata_sample_id_column for proper merging if providing a metadata file
    
    Example expression matrix format:
    ```
    transcript_id    gene_id    sample1    sample2    sample3
    ENST0001        ENSG001    100        200        150
    ENST0002        ENSG001    50         75         60
    ENST0003        ENSG002    300        250        275
    ```

    **Required Columns in Expression Matrix:**
    - `transcript_id_column_name` (default `"transcript_id"`): Identifier for each transcript
    - If `gene_id_column_name` is provided and not None, it must be a column in the expression matrix
    - All other columns are assumed to be sample expression values and must:
        - Contain numeric values
        - Have names that match the values in the metadata_sample_id_column if metadata is provided

    **Metadata Format Requirements (if provided):**
    - Must contain the metadata_sample_id_column (default "sample_id")
    - Values in this column must exactly match the sample column names in the expression matrix
    - Can contain any number of additional metadata columns

    Example metadata format:
    ```
    sample_id    condition    batch
    sample1      control      1
    sample2      treated      1
    sample3      control      2
    ```


    **Supported File Formats:**
    - `.csv`, `.tsv`, `.txt`, `.parquet`, `.xlsx` for both expression matrix and metadata files.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression matrix file. Supported file formats include `.csv`, `.tsv`, `.txt`, `.parquet`, and `.xlsx`.
    metadata_path : str, optional
        Path to the metadata file. If provided, the metadata will be merged with the expression data on the specified sample identifier column.
        Supported file formats are the same as for `expression_matrix_path`. Default is None.
    expression_measure_name : str, optional
        The name to assign to the expression measure column after melting. This will be the name of the column containing the expression values in the long-format DataFrame. Default is `"counts"`.
    cpm_normalization : bool, optional
        Whether to perform Counts Per Million (CPM) normalization on the expression data. If True, CPM values will be calculated for each sample. Default is False.
    relative_abundance : bool, optional
        Whether to calculate relative transcript abundance based on gene counts. Requires `gene_id_column_name` to be provided and not None. Default is False.
    gene_id_column_name : str, optional
        The name of the column in the expression DataFrame that contains gene identifiers. This column will remain fixed during data transformation.
        If provided and `relative_abundance` is True, relative transcript abundance will be calculated. Default is `"gene_id"`. If set to None, the gene identifier will not be used.
    transcript_id_column_name : str
        The name of the column in the expression DataFrame that contains transcript identifiers. This parameter is required and cannot be None. Default is `"transcript_id"`.
    metadata_sample_id_column : str, optional
        Column name in the metadata DataFrame that identifies samples. This column is used to merge the metadata and expression data. Default is `"sample_id"`.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame in long format containing the expression data, and optionally CPM values, relative abundances, and metadata.

    Raises
    ------
    ValueError
        If `transcript_id_column_name` is None.
        If required feature ID columns are missing in the expression matrix.
        If the expression columns are not numeric.
        If required columns are missing in the metadata DataFrame.
        If there are no overlapping sample IDs between expression data and metadata.
        If the file format is unsupported or the file cannot be read.

    Warns
    -----
    UserWarning
        If `relative_abundance` is True but `gene_id_column_name` is None.
        If there is partial overlap of sample IDs between expression data and metadata.

    Examples
    --------
    Load an expression matrix, perform CPM normalization, calculate relative transcript abundance, and merge with metadata:

    >>> from RNApysoforms import read_expression_matrix
    >>> df = read_expression_matrix(
    ...     expression_matrix_path="counts.csv",
    ...     metadata_path="metadata.csv",
    ...     expression_measure_name="counts",
    ...     cpm_normalization=True,
    ...     relative_abundance=True
    ... )
    >>> print(df.head())

    Notes
    -----
    - The `transcript_id_column_name` is set to "transcript_id" by default. The parameter is required and cannot be None.
    - The function supports multiple file formats (`.csv`, `.tsv`, `.txt`, `.parquet`, `.xlsx`) for both expression and metadata files.
    - If CPM normalization is performed, the expression measures will be scaled to reflect Counts Per Million for each sample.
    - If `gene_id_column_name` is provided and `relative_abundance` is True, relative transcript abundance is calculated as 
      `(transcript_expression / total_gene_expression) * 100`. If the total gene counts are zero, the relative abundance is set to zero to avoid division by zero errors.
    - Warnings are raised if there is only partial sample overlap between expression data and metadata.
    - The resulting DataFrame is returned in long format, with expression measures, CPM values, relative abundance for each sample-feature combination.
    - The `expression_measure_name` allows customization of the name of the expression values column in the long-format DataFrame.
    - If a metadata file is passed, the function expects that the values in the `metadata_sample_id_column` from the metadata file will be found as 
      column names in the counts matrix file.
    - Beware of using the `cpm_normalization` and `relative_abundance` options set to `True` when working with a non-raw (i.e., normalized) counts
      matrix as those results may not be accurate causing misinterpretation.
    - An example counts matrix file can be found here: https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/blob/main/tests/test_data/counts_matrix_chr21_and_Y.tsv
    - An example metadata file can be found here: https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/blob/main/tests/test_data/sample_metadata.tsv
    """

    # Check if transcript_id_column_name is None and raise an error if so
    if transcript_id_column_name is None:
        raise ValueError("The 'transcript_id_column_name' is required and cannot be None.")

    # Load the expression matrix file using the helper function
    expression_df = _get_open_file(expression_matrix_path)

    # Build a list of feature ID columns to keep fixed during data transformation
    feature_id_columns = [transcript_id_column_name]
    if gene_id_column_name is not None:
        feature_id_columns.append(gene_id_column_name)

    # Check if required feature ID columns are present in the expression DataFrame
    missing_columns = [col for col in feature_id_columns if col not in expression_df.columns]
    if missing_columns:
        raise ValueError(f"The following feature ID columns are missing in the expression dataframe: {missing_columns}")

    # Determine the expression columns by excluding feature ID columns
    expression_columns = [col for col in expression_df.columns if col not in feature_id_columns]

    # Check that expression columns are numeric
    non_numeric_columns = [col for col in expression_columns if not expression_df[col].dtype.is_numeric()]
    if non_numeric_columns:
        raise ValueError(f"The following columns are expected to be numerical but are not: {non_numeric_columns}")

    # Calculate relative transcript abundance if gene_id_column_name is provided and relative_abundance is True
    if relative_abundance and gene_id_column_name is not None:
        # Calculate total gene counts for each gene across transcripts
        expression_df = expression_df.with_columns([
            (
                pl.when(pl.col(col).sum().over(gene_id_column_name) == 0)
                .then(0)
                .otherwise((pl.col(col) / pl.col(col).sum().over(gene_id_column_name)) * 100)
                .alias(col + "_relative_abundance")
            )
            for col in expression_columns
        ])
    # If relative_abundance is True but gene_id_column_name is None, issue a warning
    elif relative_abundance and gene_id_column_name is None:
        warnings.warn(
            "relative_abundance was set to True, but gene_id_column_name was not provided (set to None). "
            "Therefore, relative abundance calculation is being skipped.",
            UserWarning
        )

    # Perform CPM normalization if requested
    if cpm_normalization:
        # Calculate total counts for each sample and scale counts to Counts Per Million (CPM)
        expression_df = expression_df.with_columns([
            (
                (pl.col(col) / pl.col(col).sum()) * 1e6
            ).alias(col + "_CPM")
            for col in expression_columns
        ])

    # Transform the expression DataFrame into long format for expression_measure_name
    expression_long = expression_df.unpivot(
        index=feature_id_columns,
        on=expression_columns,
        variable_name=metadata_sample_id_column,
        value_name=expression_measure_name
    )

    # Initialize long_expression_df with expression_long
    long_expression_df = expression_long

    # If CPM normalization was performed, melt CPM columns and join them
    if cpm_normalization:
        cpm_columns = [col + "_CPM" for col in expression_columns]

        cpm_long = expression_df.unpivot(
            index=feature_id_columns,
            on=cpm_columns,
            variable_name=metadata_sample_id_column,
            value_name="CPM"
        ).with_columns(
            # Remove the '_CPM' suffix from sample IDs if present
            pl.col(metadata_sample_id_column).str.replace(r"_CPM$", "")
        )

        # Join the CPM values to the long_expression_df
        long_expression_df = long_expression_df.join(
            cpm_long,
            on=feature_id_columns + [metadata_sample_id_column],
            how="left"
        )

    # If relative abundance was calculated, melt and join the relative abundance columns
    if relative_abundance and gene_id_column_name is not None:
        relative_abundance_columns = [col + "_relative_abundance" for col in expression_columns]

        relative_abundance_long = expression_df.unpivot(
            index=feature_id_columns,
            on=relative_abundance_columns,
            variable_name=metadata_sample_id_column,
            value_name="relative_abundance"
        ).with_columns(
            # Remove the '_relative_abundance' suffix from sample IDs if present
            pl.col(metadata_sample_id_column).str.replace(r"_relative_abundance$", "")
        )

        # Join the relative abundance values to the long_expression_df
        long_expression_df = long_expression_df.join(
            relative_abundance_long,
            on=feature_id_columns + [metadata_sample_id_column],
            how="left"
        )

    # If metadata_path is provided, load and merge metadata
    if metadata_path is not None:
        # Load the metadata file using the helper function
        metadata_df = _get_open_file(metadata_path)

        # Check if metadata_sample_id_column is present in metadata_df
        if metadata_sample_id_column not in metadata_df.columns:
            raise ValueError(
                f"The metadata_sample_id_column '{metadata_sample_id_column}' is not present in the metadata dataframe."
            )

        # Get unique sample IDs from expression data and metadata
        expression_sample_ids = long_expression_df[metadata_sample_id_column].unique().to_list()
        metadata_sample_ids = metadata_df[metadata_sample_id_column].unique().to_list()

        # Find overlapping sample IDs between expression data and metadata
        overlapping_sample_ids = set(expression_sample_ids).intersection(set(metadata_sample_ids))

        if not overlapping_sample_ids:
            raise ValueError("No overlapping sample IDs found between expression data and metadata.")

        # Warn about sample ID mismatches
        metadata_sample_ids_not_in_expression = set(metadata_sample_ids) - set(expression_sample_ids)
        expression_sample_ids_not_in_metadata = set(expression_sample_ids) - set(metadata_sample_ids)

        warning_message = ""
        if metadata_sample_ids_not_in_expression:
            warning_message += (f"The following sample IDs are present in metadata but not in expression data: {list(metadata_sample_ids_not_in_expression)}. "
                                "Only returning sample IDs that are present in both.")
        if expression_sample_ids_not_in_metadata:
            warning_message += (f"The following sample IDs are present in expression data but not in metadata: {list(expression_sample_ids_not_in_metadata)}. "
                                "Only returning sample IDs that are present in both.")
        if warning_message:
            warnings.warn(warning_message)

        # Merge metadata with the long_expression_df on metadata_sample_id_column
        long_expression_df = long_expression_df.join(
            metadata_df,
            on=metadata_sample_id_column,
            how="inner"
        )

    # Return the final long-format DataFrame
    return long_expression_df

def _get_open_file(file_path: str) -> pl.DataFrame:
    """
    Opens a file based on its extension and loads it into a Polars DataFrame.

    This helper function supports multiple file formats such as `.csv`, `.tsv`, `.txt`, `.parquet`, and `.xlsx`.
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
        If the file extension is unsupported.
        If the file cannot be read due to an error.

    Examples
    --------
    Open a CSV file:

    >>> df = _get_open_file("data.csv")

    Open a TSV file:

    >>> df = _get_open_file("data.tsv")

    Notes
    -----
    - This function is used internally by `read_expression_matrix` to load expression and metadata files.
    - It handles different file extensions and raises a clear error if the file format is unsupported.
    """

    # Extract the file extension to determine the file format
    _, file_extension = os.path.splitext(file_path)

    try:
        # Open the file based on its extension
        if file_extension in [".tsv", ".txt"]:
            # Read tab-separated values
            return pl.read_csv(file_path, separator="\t")
        elif file_extension == ".csv":
            # Read comma-separated values
            return pl.read_csv(file_path)
        elif file_extension == ".parquet":
            # Read Parquet file
            return pl.read_parquet(file_path)
        elif file_extension == ".xlsx":
            # Read Excel file
            return pl.read_excel(file_path)
        else:
            # Raise an error for unsupported file extensions
            raise ValueError(
                f"Unsupported file extension '{file_extension}'. Supported extensions are .tsv, .txt, .csv, .parquet, .xlsx"
            )
    except Exception as e:
        # Raise an error if the file cannot be read
        raise ValueError(f"Failed to read the file '{file_path}': {e}")
