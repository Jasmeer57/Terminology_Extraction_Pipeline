rule generate_embeddings:
    input:
        "results/linked_entities.json"
    output:
        "results/embeddings_comparison.csv"
    params:
        models=lambda wildcards, config: config["huggingface_models"]
    script:
        "scripts/generate_embeddings.py"
