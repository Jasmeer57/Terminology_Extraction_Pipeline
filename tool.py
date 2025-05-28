import tkinter as tk
from tkinter import filedialog, messagebox
import os
import pandas as pd
import json
import requests
from Bio import Entrez
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModel
from sentence_transformers import SentenceTransformer
from sklearn.cluster import DBSCAN

# === Configuration ===
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# HuggingFace Models for Comparison
HUGGINGFACE_MODELS = {
    "BioBERT": "dmis-lab/biobert-base-cased-v1.2",
    "PubMedBERT": "microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext",
    "SciBERT": "allenai/scibert_scivocab_cased",
    "BioClinicalBERT": "emilyalsentzer/Bio_ClinicalBERT",
    "SapBERT": "cambridgeltl/sapbert-from-pubmed",
    "BioLinkBERT": "microsoft/BiomedNLP-LinkBERT-base-4096"
}

# === Pipeline Functions ===
def search_pubmed(query, email, max_results=1000):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_abstracts(pmid_list, email, batch_size=200):
    Entrez.email = email
    abstracts = []
    for start in tqdm(range(0, len(pmid_list), batch_size)):
        end = min(start + batch_size, len(pmid_list))
        batch_pmids = pmid_list[start:end]
        ids = ",".join(batch_pmids)

        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for article in records["PubmedArticle"]:
            pmid = article["MedlineCitation"]["PMID"]
            article_data = article["MedlineCitation"]["Article"]
            title = article_data.get("ArticleTitle", "")
            abstract_text = ""
            if "Abstract" in article_data:
                abstract_obj = article_data["Abstract"]
                abstract_items = abstract_obj.get("AbstractText", [])
                if isinstance(abstract_items, list):
                    abstract_text = " ".join(str(item) for item in abstract_items)
                else:
                    abstract_text = str(abstract_items)
            abstracts.append({"pmid": pmid, "title": title, "abstract": abstract_text})
    return abstracts

def extract_terms(abstracts):
    terms = set()
    for abstract in abstracts:
        words = abstract.split()
        terms.update(words)
    return list(terms)

def entity_linking(term):
    try:
        res = requests.get(f"https://www.ebi.ac.uk/ols/api/search?q={term}&ontology=envo")
        if res.status_code == 200:
            hits = res.json()["response"]["docs"]
            return hits[0]["label"] if hits else None
    except Exception:
        return None

def compare_models(terms):
    results = {}
    for model_name, model_id in HUGGINGFACE_MODELS.items():
        try:
            tokenizer = AutoTokenizer.from_pretrained(model_id)
            model = AutoModel.from_pretrained(model_id)
            embeddings = [tokenizer(term, return_tensors="pt") for term in terms]
            results[model_name] = embeddings
        except Exception as e:
            print(f"Error with {model_name}: {e}")
    return results

def cluster_terms(terms):
    model = SentenceTransformer("all-MiniLM-L6-v2")
    embeddings = model.encode(terms, show_progress_bar=True)
    clustering = DBSCAN(eps=1.5, min_samples=3, metric="cosine").fit(embeddings)
    return pd.DataFrame({"term": terms, "cluster": clustering.labels_})

def run_pipeline(search_term, email):
    try:
        print("🔍 Searching PubMed...")
        pmids = search_pubmed(search_term, email)
        print(f"Found {len(pmids)} PubMed articles.")

        print("📄 Fetching abstracts...")
        abstracts = fetch_abstracts(pmids, email)
        df_abstracts = pd.DataFrame(abstracts)
        df_abstracts.to_csv(os.path.join(OUTPUT_DIR, "pubmed_abstracts.csv"), index=False)

        print("🔎 Extracting terms...")
        terms = extract_terms(df_abstracts["abstract"].tolist())

        print("🔗 Performing entity linking...")
        linked_terms = [{"term": term, "linked_term": entity_linking(term)} for term in terms]
        df_linked_terms = pd.DataFrame(linked_terms)
        df_linked_terms.to_csv(os.path.join(OUTPUT_DIR, "linked_terms.csv"), index=False)

        print("📊 Comparing models...")
        model_results = compare_models(terms)
        with open(os.path.join(OUTPUT_DIR, "model_comparison.json"), "w") as f:
            json.dump(model_results, f, indent=2)

        print("📈 Clustering terms...")
        clustered_terms = cluster_terms(terms)
        clustered_terms.to_csv(os.path.join(OUTPUT_DIR, "clustered_terms.csv"), index=False)

        messagebox.showinfo("Pipeline Complete", f"Pipeline executed successfully!\nResults saved to {OUTPUT_DIR}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

# === GUI Setup ===
def start_pipeline():
    search_term = search_term_var.get().strip()
    email = email_var.get().strip()
    if not search_term:
        messagebox.showwarning("Missing Input", "Please enter a search string.")
        return
    if not email:
        messagebox.showwarning("Missing Input", "Please enter your email ID.")
        return
    run_pipeline(search_term, email)

app = tk.Tk()
app.title("PubMed Pipeline Tool")
app.geometry("500x350")

# Input Search String
tk.Label(app, text="Enter PubMed Search String:").pack(pady=5)
search_term_var = tk.StringVar()
tk.Entry(app, textvariable=search_term_var, width=50).pack(pady=5)

# Input Email ID
tk.Label(app, text="Enter Your Email ID:").pack(pady=5)
email_var = tk.StringVar()
tk.Entry(app, textvariable=email_var, width=50).pack(pady=5)

# Start Pipeline Button
tk.Button(app, text="Run Pipeline", command=start_pipeline, bg="green", fg="white").pack(pady=20)

# Run the application
app.mainloop()