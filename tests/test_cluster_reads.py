import pandas as pd
from src.pipeline import cluster_consolidated_reads_by_edit_distance

def test_basic_clustering():
    df = pd.DataFrame({
        "read_sequence": ["GGAAGGAAGGAAGGAA", "GGGAGGAA", "GGAAGGAAGGAAGGAT", "GGAAGGACGGAT"],
        "read_support": [8, 4, 2, 2]
    })

    # Use a clustering threshold of 4 (only edit distance ≤ 4 allowed)
    clustered = cluster_consolidated_reads_by_edit_distance(df, cluster_dist=4)

    # confirm that the cluster id column was generated from the function
    assert "cluster_id" in clustered.columns
    assert "read_sequence" in clustered.columns

    # Check number of clusters
    assert len(clustered) == 3

    # Ensure read_support was added when reads were comvined into a single cluster
    top_cluster = clustered.sort_values("read_support", ascending=False).iloc[0]
    assert top_cluster["read_support"] == 10  # 8 + 2 from similar reads

