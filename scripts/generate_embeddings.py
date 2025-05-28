import json
import pandas as pd
from sentence_transformers import SentenceTransformer, util
from snakemake import input, output, params

with open(input[0], "r") as f:
    linked_terms = json.load(f)

terms = list(linked_terms.keys())
model_names = params.models
model_scores = []

# Load models and compute embeddings
for model_name in model_names:
    model = SentenceTransformer(model_name)
    embeddings = model.encode(terms, convert_to_tensor=True)
    
    # Pairwise cosine similarity (average as a proxy for internal consistency)
    similarity = util.pytorch_cos_sim(embeddings, embeddings)
    avg_sim = similarity.mean().item()
    
    model_scores.append({
        "model": model_name,
        "average_similarity": round(avg_sim, 4)
    })

df = pd.DataFrame(model_scores)
df.to_csv(output[0], index=False)
