import polars as pl
import warnings
from typing import Union
from RNApysoforms.utils import check_df

def gene_filtering(
    target_gene: str,
    annotation: pl.DataFrame,
    expression_matrix: pl.DataFrame = None,
    transcript_id_column: str = "transcript_id",
    gene_id_column: str = "gene_name",
    order_by_expression_column: str = "counts",
    order_by_expression: bool = True,
    keep_top_expressed_transcripts: Union[str, int] = "all"
) -> Union[pl.DataFrame, tuple]:
    """
    Filters genomic annotations and optionally an expression matrix for a specific gene,
    with options to order and select top expressed transcripts.

    This function filters the provided annotation DataFrame to include only entries corresponding
    to the specified `target_gene`, identified using the column specified by `gene_id_column`.
    If an expression matrix is provided, it will also be filtered to retain only the entries
    corresponding to the filtered transcripts based on the `transcript_id_column`. Additionally,
    it provides options to order transcripts by their total expression levels and to keep only
    the top expressed transcripts, specified by `keep_top_expressed_transcripts`.

    **Required Columns in `annotation` DataFrame:**
    - `gene_id_column` (default `"gene_name"`): Column containing gene identifiers used for filtering.
    - `transcript_id_column` (default `"transcript_id"`): Column containing transcript identifiers.

    **Required Columns in `expression_matrix` DataFrame (if provided):**
    - `transcript_id_column` (same as in `annotation`): Column containing transcript identifiers matching those in `annotation`.
    - `order_by_expression_column` (default `"counts"`): Column containing expression values used for ordering and filtering.

    Parameters
    ----------
    target_gene : str
        The gene identifier to filter in the annotation DataFrame.
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations. Must include the columns specified by `gene_id_column`
        and `transcript_id_column`.
    expression_matrix : pl.DataFrame, optional
        A Polars DataFrame containing expression data. If provided, it will be filtered to match the filtered
        annotation based on `transcript_id_column`. Default is None.
    transcript_id_column : str, optional
        The column name representing transcript identifiers in both the annotation and expression matrix.
        Default is 'transcript_id'.
    gene_id_column : str, optional
        The column name in the annotation DataFrame that contains gene identifiers used for filtering.
        Default is 'gene_name'.
    order_by_expression_column : str, optional
        The column name in the expression matrix that contains expression values used for ordering and filtering.
        Default is 'counts'.
    order_by_expression : bool, optional
        If True, transcripts will be ordered by their total expression levels in descending order.
        Default is True.
    keep_top_expressed_transcripts : Union[str, int], optional
        Determines the number of top expressed transcripts to keep after ordering by expression levels.
        Can be 'all' to keep all transcripts or an integer to keep the top N transcripts. Default is 'all'.

    Returns
    -------
    pl.DataFrame or tuple
        - If `expression_matrix` is provided, returns a tuple of (filtered_annotation, filtered_expression_matrix).
        - If `expression_matrix` is None, returns only the `filtered_annotation`.

    Raises
    ------
    TypeError
        If `annotation` or `expression_matrix` are not Polars DataFrames.
    ValueError
        If required columns are missing in the `annotation` or `expression_matrix` DataFrames.
    ValueError
        If the filtered expression matrix is empty after filtering.
    ValueError
        If `keep_top_expressed_transcripts` is not 'all' or a positive integer.
    Warning
        If there are transcripts present in the annotation but missing in the expression matrix.

    Examples
    --------
    Filter an annotation DataFrame by a specific gene:

    >>> import polars as pl
    >>> from RNApysoforms.annotation import gene_filtering
    >>> annotation_df = pl.DataFrame({
    ...    "gene_name": ["APP", "APP", "APP"],
    ...    "transcript_id": ["tx1", "tx2", "tx3"]
    ... })
    >>> expression_matrix_df = pl.DataFrame({
    ...    "transcript_id": ["tx1", "tx2", "tx3"],
    ...    "counts": [300, 100, 200]
    ... })
    >>> target_gene = "APP"
    >>> filtered_annotation, filtered_expression_matrix = gene_filtering(
    ...    target_gene,
    ...    annotation_df,
    ...    expression_matrix=expression_matrix_df,
    ...    order_by_expression=True
    ... )

    Notes
    -----
    - The function filters the `annotation` DataFrame to include only entries where `gene_id_column` matches `target_gene`.
    - If an `expression_matrix` is provided, the function filters it to include only transcripts present in the filtered annotation.
    - The function checks for transcripts present in the annotation but missing in the expression matrix and issues a warning for such discrepancies.
    - If `order_by_expression` is True, transcripts are ordered by their total expression levels computed from the `order_by_expression_column` in the expression matrix.
    - If `keep_top_expressed_transcripts` is an integer, only the top N expressed transcripts are kept after ordering.
    - If `keep_top_expressed_transcripts` is 'all', all transcripts are kept.
    - If transcripts are present in the expression matrix but not in the annotation, they are silently ignored, and only overlapping transcripts are returned without a warning.

    """

    # Validate the input 'annotation' DataFrame
    # Check if 'annotation' is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}."
            "\nYou can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
        )

    # Ensure required columns are present in the 'annotation' DataFrame
    check_df(annotation, [gene_id_column, transcript_id_column])

    # Filter the annotation DataFrame to include only entries for the target gene
    filtered_annotation = annotation.filter(pl.col(gene_id_column) == target_gene)

    # If no entries are found for the target gene, raise a ValueError
    if filtered_annotation.is_empty():
        raise ValueError(f"No annotation found for gene: {target_gene} in the '{gene_id_column}' column")

    if expression_matrix is not None:
        # Validate the input 'expression_matrix' DataFrame
        # Check if 'expression_matrix' is a Polars DataFrame
        if not isinstance(expression_matrix, pl.DataFrame):
            raise TypeError(
                f"Expected 'expression_matrix' to be of type pl.DataFrame, got {type(expression_matrix)}."
                "\nYou can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
            )

        # Ensure required columns are present in the expression matrix
        check_df(expression_matrix, [transcript_id_column, order_by_expression_column])

        # Filter the expression matrix to include only transcripts present in the filtered annotation
        filtered_expression_matrix = expression_matrix.filter(
            pl.col(transcript_id_column).is_in(filtered_annotation[transcript_id_column].to_list())
        )

        # If the filtered expression matrix is empty after filtering, raise a ValueError
        if filtered_expression_matrix.is_empty():
            raise ValueError(
                f"Expression matrix is empty after filtering. No matching '{transcript_id_column}' entries "
                f"between expression matrix and annotation found for gene '{target_gene}'."
            )

        # Identify transcripts present in annotation but missing in expression matrix
        # Get sets of transcripts in annotation and expression matrix
        annotation_transcripts = set(filtered_annotation[transcript_id_column].unique())
        expression_transcripts = set(filtered_expression_matrix[transcript_id_column].unique())

        # Find transcripts that are in annotation but not in expression matrix
        missing_in_expression = annotation_transcripts - expression_transcripts

        # Transcripts present in expression matrix but not in annotation are silently ignored

        # Issue a warning for transcripts missing in the expression matrix
        if missing_in_expression:
            warnings.warn(
                f"{len(missing_in_expression)} transcript(s) are present in the annotation but missing in the expression matrix. "
                f"Missing transcripts: {', '.join(sorted(missing_in_expression))}. "
                "Only transcripts present in both will be returned."
            )

        # Ensure both filtered_annotation and filtered_expression_matrix contain only common transcripts
        common_transcripts = annotation_transcripts & expression_transcripts
        filtered_annotation = filtered_annotation.filter(
            pl.col(transcript_id_column).is_in(list(common_transcripts))
        )
        filtered_expression_matrix = filtered_expression_matrix.filter(
            pl.col(transcript_id_column).is_in(list(common_transcripts))
        )

        # Aggregate expression data to compute total expression per transcript
        aggregated_df = filtered_expression_matrix.group_by(transcript_id_column).agg(
            pl.col(order_by_expression_column).sum().alias("total_expression")
        )

        # Sort transcripts by total expression in descending order
        sorted_transcripts = aggregated_df.sort("total_expression", descending=True)

        if order_by_expression:

            # Order the filtered_annotation and filtered_expression_matrix by total expression
            # Join total expression back to annotation and expression matrix
            filtered_annotation = filtered_annotation.join(
                sorted_transcripts.select([transcript_id_column, "total_expression"]),
                on=transcript_id_column,
                how="inner"
            ).sort("total_expression", descending=False).drop("total_expression")

            filtered_expression_matrix = filtered_expression_matrix.join(
                sorted_transcripts.select([transcript_id_column, "total_expression"]),
                on=transcript_id_column,
                how="inner"
            ).sort("total_expression", descending=False).drop("total_expression")

        # Determine transcripts to keep based on 'keep_top_expressed_transcripts'
        if isinstance(keep_top_expressed_transcripts, int) and keep_top_expressed_transcripts > 0:
            # Keep only the top N expressed transcripts
            if keep_top_expressed_transcripts < len(sorted_transcripts):
                transcripts_to_keep = sorted_transcripts.head(keep_top_expressed_transcripts)[transcript_id_column]
            else:
                # If requested number exceeds available transcripts, keep all and issue a warning
                transcripts_to_keep = sorted_transcripts[transcript_id_column]
                warnings.warn(
                    "The number specified in 'keep_top_expressed_transcripts' exceeds the total number of transcripts. "
                    "All transcripts will be kept."
                )
        elif keep_top_expressed_transcripts == "all":
            # Keep all transcripts
            transcripts_to_keep = sorted_transcripts[transcript_id_column]
        else:
            # Raise error if 'keep_top_expressed_transcripts' is invalid
            raise ValueError(
                f"'keep_top_expressed_transcripts' must be 'all' or a positive integer, got {keep_top_expressed_transcripts}."
            )

        # Filter annotation and expression matrix to include only the selected transcripts
        filtered_annotation = filtered_annotation.filter(
            pl.col(transcript_id_column).is_in(transcripts_to_keep.to_list())
        )
        filtered_expression_matrix = filtered_expression_matrix.filter(
            pl.col(transcript_id_column).is_in(transcripts_to_keep.to_list())
        )

        # Return the filtered annotation and expression matrix
        return filtered_annotation, filtered_expression_matrix

    else:
        # If no expression_matrix is provided, return only the filtered annotation
        return filtered_annotation
