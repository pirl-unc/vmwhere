import pandas as pd
from vmwhere.genotyper import cluster_consolidated_reads_by_edit_distance

def test_basic_clustering():
    df = pd.DataFrame({
        "ALT": ["GGAAGGAAGGAAGGAA", "GGGAGGAA", "GGAAGGAAGGAAGGAT", "GGAAGGACGGAT"],
        "RS": [8, 4, 2, 2]
    })

    # Use a clustering threshold of 4 (only edit distance ≤ 4 allowed)
    clustered = cluster_consolidated_reads_by_edit_distance(df, cluster_dist=4)

    # confirm that the expected columns exist
    assert "ALT" in clustered.columns
    assert "RS" in clustered.columns

    # Check number of clusters
    assert len(clustered) == 3

    # Ensure RS was summed when reads were combined into a single cluster
    top_cluster = clustered.sort_values("RS", ascending=False).iloc[0]
    assert top_cluster["RS"] == 10  # 8 + 2 from similar reads

