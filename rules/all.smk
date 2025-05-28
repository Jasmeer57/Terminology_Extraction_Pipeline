rule all:
    input:
        "results/abstracts.json",
        "results/extracted_terms.json",
        "results/linked_entities.json",
        "results/embeddings_comparison.csv",
        "results/clusters.json"
