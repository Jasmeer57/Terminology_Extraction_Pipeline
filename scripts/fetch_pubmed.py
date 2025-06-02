import argparse
from Bio import Entrez
import json
import os

def fetch_pubmed(email, search_term, output_file):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=10)
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]
    if not ids:
        print("No PubMed IDs found for search term:", search_term)
        return

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    abstracts = []
    for article in records["PubmedArticle"]:
        medline = article["MedlineCitation"]
        pmid = medline["PMID"]
        article_data = medline["Article"]
        title = article_data.get("ArticleTitle", "")
        abstract_text = ""
        if "Abstract" in article_data and "AbstractText" in article_data["Abstract"]:
            abstract_text = " ".join(article_data["Abstract"]["AbstractText"])
        abstracts.append({
            "PMID": str(pmid),
            "Title": title,
            "Abstract": abstract_text
        })

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(abstracts, f, indent=2)

    print(f"Saved {len(abstracts)} abstracts to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubMed abstracts")
    parser.add_argument("--email", required=True, help="Email address for NCBI Entrez")
    parser.add_argument("--search_term", required=True, help="Search term for PubMed query")
    parser.add_argument("--output", default="results/abstracts.json", help="Output JSON file")
    args = parser.parse_args()

    fetch_pubmed(args.email, args.search_term, args.output)
