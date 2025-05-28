import json
import numpy as np
from sentence_transformers import SentenceTransformer
from sklearn.cluster import DBSCAN
from snakemake import input, output, config

# Load terms and model
with open("results/linked_entities.json", "r") as f:
    linked_terms = json.load(f)

terms = list(linked_terms.keys())
model_name = config["huggingface_models"][0]
model = SentenceTransformer(model_name)

embeddings = model.encode(terms)
clustering = DBSCAN(eps=1.0, min_samples=2).fit(embeddings)

clusters = {}
for term, label in zip(terms, clustering.labels_):
    clusters.setdefault(str(label), []).append(term)

with open(output[0], "w") as f:
    json.dump(clusters, f, indent=2)
