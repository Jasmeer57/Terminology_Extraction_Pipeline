rule all:
    input:
        "results/abstracts.json",
        "results/extracted_terms.json",
        "results/linked_entities.json",
        "results/embeddings_comparison.csv",
        "results/clusters.json"
        
print("Email:", config["email"])
print("Search term:", config["search_term"])
print("Output dir:", config["output_dir"])
print("Models:", config["huggingface_models"])
