rule embeddings:
    input:
        "results/linked_entities.json"
    output:
        "results/embeddings_comparison.csv"
    params:
        models=" ".join(config["huggingface_models"])
    shell:
        """
        python scripts/compare_embeddings.py --input {input} --output {output} --models {params.models}
        """
