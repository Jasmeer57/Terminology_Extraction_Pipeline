rule cluster_terms:
    input:
        "results/embeddings_comparison.csv"
    output:
        "results/clusters.json"
    script:
        "scripts/cluster_terms.py"
