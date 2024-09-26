import pandas as pd
from typing import Union, List

def to_intron(exons: pd.DataFrame, group_var: Union[str, List[str]] = None) -> pd.DataFrame:
    """
    Convert exon coordinates to intron coordinates.

    Parameters:
    -----------
    exons : pd.DataFrame
        DataFrame containing exon coordinates. Must have 'start' and 'end' columns.
    group_var : str or List[str], optional
        Column name(s) for grouping transcripts. Default is None.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the intron coordinates.
    """
    # Input checks
    required_cols = ['start', 'end']
    assert all(col in exons.columns for col in required_cols), \
        f"exons DataFrame must have columns: {', '.join(required_cols)}"
    
    if group_var:
        if isinstance(group_var, str):
            group_var = [group_var]
        assert all(col in exons.columns for col in group_var), \
            f"group_var columns must exist in exons DataFrame"

    # Sort exons by start and end coordinates
    sort_cols = group_var + ['start', 'end'] if group_var else ['start', 'end']
    exons_sorted = exons.sort_values(sort_cols)

    # Group by transcript if group_var is provided
    if group_var:
        grouped = exons_sorted.groupby(group_var)
    else:
        # If no group_var, treat all exons as one group
        grouped = [(None, exons_sorted)]

    introns_list = []

    for _, group in grouped:
        # Calculate intron start and end
        intron_start = group['end'].shift()
        intron_end = group['start']

        # Create introns DataFrame
        introns = pd.DataFrame({
            'start': intron_start,
            'end': intron_end,
            'type': 'intron'
        })

        # Add group variables if present
        if group_var:
            for var in group_var:
                introns[var] = group[var].iloc[0]

        # Add other columns from exons, if present
        for col in exons.columns:
            if col not in introns.columns and col not in ['start', 'end']:
                introns[col] = group[col].iloc[0]

        introns_list.append(introns)

    # Combine all introns
    all_introns = pd.concat(introns_list, ignore_index=True)

    # Remove NAs and introns with width of 1
    all_introns = all_introns.dropna(subset=['start', 'end'])
    all_introns = all_introns[abs(all_introns['end'] - all_introns['start']) != 1]

    return all_introns

# Example usage
#if __name__ == "__main__":
    # Create sample exon data
#    exons = pd.DataFrame({
#        'start': [100, 300, 600, 150, 400, 700],
#        'end': [200, 400, 700, 250, 500, 800],
#        'strand': ['+', '+', '+', '-', '-', '-'],
#        'seqnames': ['chr1'] * 6,
#        'transcript_name': ['tx1', 'tx1', 'tx1', 'tx2', 'tx2', 'tx2'],
#        'type': ['exon'] * 6
#    })

    # Convert exons to introns
 #   introns = to_intron(exons, group_var='transcript_name')

 #   print("Exons:")
 #   print(exons)
 #   print("\nIntrons:")
 #   print(introns)