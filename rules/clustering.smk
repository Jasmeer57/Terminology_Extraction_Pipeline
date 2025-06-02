rule clustering:
    input:
        "results/embeddings_comparison.csv"
    output:
        "results/clusters.json"
    shell:
        """
        python scripts/cluster_terms.py --input {input} --output {output}
        """
