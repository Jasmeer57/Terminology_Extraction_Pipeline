
# Snakefile

# Read user config from command line
email = config["email"]
search_term = config["search_term"]

# Set default config variables
config["output_dir"] = "results/"
config["huggingface_models"] = [
    "dmis-lab/biobert-base-cased-v1.1",
    "allenai/scibert_scivocab_uncased"
]

# Debug prints to check config values
print("Email:", email)
print("Search term:", search_term)
print("Output dir:", config["output_dir"])
print("Models:", config["huggingface_models"])

# Include all rule files
include: "rules/all.smk"
include: "rules/fetch_pubmed.smk"
include: "rules/extract_terms.smk"
include: "rules/entity_linking.smk"
include: "rules/embeddings.smk"
include: "rules/clustering.smk"
