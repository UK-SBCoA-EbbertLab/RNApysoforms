import polars as pl
import warnings
from RNApysoforms.utils import check_df

def gene_filtering(
    target_gene: str,
    annotation: pl.DataFrame,
    expression_matrix: pl.DataFrame = None,
    group_var: str = "transcript_id",
    gene_id_column: str = "gene_name"
):
    """
    Filters genomic annotations and, optionally, an expression matrix based on a specific gene identifier.

    This function filters the provided annotation DataFrame to include only the specified gene, identified by the `target_gene`.
    The gene is identified using the column specified by `gene_id_column`. If an expression matrix is provided, it will also
    filter the expression matrix to retain only the entries corresponding to the filtered gene based on the `group_var` column.
    The function checks for missing entries between the annotation and expression matrix and issues warnings if discrepancies
    are found.

    Parameters
    ----------
    target_gene : str
        The gene identifier to filter in the annotation DataFrame.
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations. Must include the column specified by `gene_id_column` and the column specified by `group_var`.
    expression_matrix : pl.DataFrame, optional
        A Polars DataFrame containing expression data. If provided, it will be filtered to match the filtered annotation based on `group_var`.
        Default is None.
    group_var : str, optional
        The column name representing group identifiers (e.g., transcript IDs) in both the annotation and expression matrix.
        Default is 'transcript_id'.
    gene_id_column : str, optional
        The column name in the annotation DataFrame that contains gene identifiers used for filtering.
        Default is 'gene_name'.

    Returns
    -------
    tuple or pl.DataFrame
        - If `expression_matrix` is provided, returns a tuple of (filtered_annotation, filtered_expression_matrix).
        - If `expression_matrix` is None, returns only the `filtered_annotation`.

    Raises
    ------
    TypeError
        If the `annotation` or `expression_matrix` are not Polars DataFrames.
    ValueError
        If required columns (`gene_id_column`, `group_var`) are missing in the annotation DataFrame,
        or if the expression matrix does not contain the `group_var`.
    ValueError
        If the filtered expression matrix is empty after filtering.
    Warning
        If there are discrepancies between groups in the annotation and expression matrix.

    Examples
    --------
    Filter an annotation DataFrame by a specific gene:

    >>> import polars as pl
    >>> annotation_df = pl.DataFrame({
    ...     "gene_name": ["BRCA1", "BRCA1", "TP53"],
    ...     "transcript_id": ["tx1", "tx2", "tx3"]
    ... })
    >>> filtered_annotation = gene_filtering("BRCA1", annotation_df)
    >>> print(filtered_annotation)

    Filter an annotation and expression matrix by a specific gene:

    >>> expression_matrix_df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx2", "tx4"],
    ...     "sample1": [100, 200, 300],
    ...     "sample2": [150, 250, 350]
    ... })
    >>> filtered_annotation, filtered_expression_matrix = gene_filtering(
    ...     "BRCA1", annotation_df, expression_matrix=expression_matrix_df
    ... )
    >>> print(filtered_annotation)
    >>> print(filtered_expression_matrix)

    Notes
    -----
    - The function assumes the `annotation` DataFrame contains the columns specified by `gene_id_column` and `group_var`.
    - If an `expression_matrix` is provided, the function checks for discrepancies between groups in the annotation and
      expression matrix and raises warnings if there are differences.
    - If no matching groups are found after filtering the expression matrix, the function raises an error.
    """

    # Check if annotation is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected annotation to be of type pl.DataFrame, got {type(annotation)}" +
            "\n You can use polars_df = pandas_df.from_pandas() to convert a pandas df into a polars df"
        )

    # Check if annotation has the specified gene_id_column and group_var columns
    check_df(annotation, [gene_id_column, group_var])

    # Filter annotation based on 'target_gene'
    filtered_annotation = annotation.filter(pl.col(gene_id_column) == target_gene)

    # Check if filtered_annotation is empty and raise an error if true
    if filtered_annotation.is_empty():
        raise ValueError(f"No annotation found for gene: {target_gene} in the '{gene_id_column}' column")

    if expression_matrix is not None:
        # Check if expression_matrix is a Polars DataFrame
        if not isinstance(expression_matrix, pl.DataFrame):
            raise TypeError("Counts matrix must be a polars DataFrame.")

        # Check if expression_matrix has the group_var column
        check_df(expression_matrix, [group_var])

        # Filter expression_matrix based on groups in filtered_annotation
        filtered_expression_matrix = expression_matrix.filter(
            pl.col(group_var).is_in(filtered_annotation[group_var])
        )

        # If filtered expression matrix is empty, throw informative error
        if filtered_expression_matrix.is_empty():
            raise ValueError(
                f"Filtered expression matrix is empty after filtering. No matching entries {group_var} entries found for {target_gene} gene."
            )

        # Check for discrepancies between groups in annotation and expression_matrix
        annotation_groups = set(filtered_annotation[group_var].unique())
        expression_groups = set(filtered_expression_matrix[group_var].unique())

        # Groups in annotation but not in expression matrix
        missing_in_expression = annotation_groups - expression_groups

        # Groups in expression matrix but not in annotation
        missing_in_annotation = expression_groups - annotation_groups

        # Warn about groups missing in the expression matrix
        if missing_in_expression:
            warnings.warn(
                f"{len(missing_in_expression)} group(s) are present in the annotation but missing in the expression matrix. "
                f"Missing groups: {', '.join(sorted(missing_in_expression))}"
            )

        # Warn about groups missing in the annotation
        if missing_in_annotation:
            warnings.warn(
                f"{len(missing_in_annotation)} group(s) are present in the expression matrix but missing in the annotation. "
                f"Missing groups: {', '.join(sorted(missing_in_annotation))}"
            )

        return filtered_annotation, filtered_expression_matrix
    else:
        return filtered_annotation
