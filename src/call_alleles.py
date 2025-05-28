

def identify_alleles(read_summary_df, minor_thresh=0.15, homozygous_thresh=0.8):
    """Determine alleles based on read support."""

    read_summary_df['total_region_reads'] = read_summary_df.groupby('region_id')['read_support'].transform('sum')
    read_summary_df['percent_support'] = read_summary_df['read_support'] / read_summary_df['total_region_reads']

    allele_dfs = []

    for region_id, group in read_summary_df.groupby('region_id'):
        max_support = group['percent_support'].max()
        if max_support >= homozygous_thresh:
            filtered = group[group['percent_support'] >= homozygous_thresh]
        else:
            filtered = group[group['percent_support'] >= minor_thresh]
        allele_dfs.append(filtered)

    return pd.concat(allele_dfs, ignore_index=True)

