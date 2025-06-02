import json
import pandas as pd
from sentence_transformers import SentenceTransformer, util
import sys

def compute_embeddings(terms, model_name):
    model = SentenceTransformer(model_name)
    embeddings = model.encode(terms, convert_to_tensor=True)
    return embeddings

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute term embeddings and similarity.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model", required=True)
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        terms = json.load(f)

    embeddings = compute_embeddings(terms, args.model)

    # Compute cosine similarity matrix
    cos_sim_matrix = util.pytorch_cos_sim(embeddings, embeddings).cpu().numpy()

    # Convert to DataFrame
    df = pd.DataFrame(cos_sim_matrix, index=terms, columns=terms)
    df.to_csv(args.output)
