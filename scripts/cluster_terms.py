import json
import pandas as pd
from sklearn.cluster import DBSCAN
import numpy as np
import sys

def cluster_embeddings(similarity_csv, eps=0.5, min_samples=2):
    df = pd.read_csv(similarity_csv, index_col=0)
    # Convert similarity to distance
    distance_matrix = 1 - df.values
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = clustering.fit_predict(distance_matrix)
    return dict(zip(df.index, labels.tolist()))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Cluster terms using DBSCAN.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    clusters = cluster_embeddings(args.input)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(clusters, f, indent=2)
