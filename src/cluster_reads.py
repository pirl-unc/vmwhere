

def cluster_consolidated_reads_by_edit_distance(consolidated_results, cluster_dist):
    """
    Clusters read_sequences based on Levenshtein distance and read support.

    Parameters:
        df (pd.DataFrame): Must contain 'read_sequence' and 'read_support'.
        cluster_distance (int): Max allowed distance to join a cluster.

    Returns:
        pd.DataFrame: DataFrame with cluster assignment and updated supporting reads for the cluster
    """
    results_sorted = consolidated_results.sort_values(by='read_support', ascending=False).reset_index(drop=True)
    clusters = []  # Each cluster is {'centroid': str, 'members': [int], 'ids': [int]}
    cluster_assignments = [-1] * len(results_sorted)

    for i, row in results_sorted.iterrows():
        sequence = row['read_sequence']
        best_dist = float('inf')
        best_cluster_idx = -1

        # Compare to each existing cluster centroid
        for idx, cluster in enumerate(clusters):
            dist = Levenshtein.distance(sequence, cluster['centroid'])
            if dist < best_dist and dist <= cluster_dist:
                best_dist = dist
                best_cluster_idx = idx

        if best_cluster_idx == -1:
            # Create new cluster with this sequence as centroid
            clusters.append({'centroid': sequence, 'members': [i]})
            cluster_assignments[i] = len(clusters) - 1
        else:
            # Assign to closest cluster
            clusters[best_cluster_idx]['members'].append(i)
            cluster_assignments[i] = best_cluster_idx

    clustered = results_sorted.copy()
    clustered['cluster_id'] = cluster_assignments
    # identify clustered total read support for each cluster
    cluster_read_support = clustered.groupby('cluster_id')['read_support'].sum()
    # identify the index of the centroid for each cluster
    centroids_idx = clustered.groupby('cluster_id')['read_support'].idxmax()
    # create dataframe with only centroids and drop read support column since it's wrong
    centroids_df = clustered.loc[centroids_idx]
    centroids_df = centroids_df.drop(columns=['read_support'])
    # update total read support for each cluster
    centroids_df_final = centroids_df.merge(cluster_read_support, on='cluster_id', how='left')

    return centroids_df_final

