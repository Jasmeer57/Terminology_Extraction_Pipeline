import requests
import json
import time
import sys

def link_terms_to_ontology(terms):
    base_url = "https://www.ebi.ac.uk/ols/api/search"
    linked = []

    for term in terms:
        params = {
            "q": term,
            "ontology": "all",
            "type": "class",
            "size": 1
        }
        r = requests.get(base_url, params=params)
        if r.status_code == 200:
            results = r.json().get('response', {}).get('docs', [])
            if results:
                linked.append({
                    "term": term,
                    "ontology_term": results[0].get('label'),
                    "ontology_iri": results[0].get('iri')
                })
        time.sleep(0.1)
    return linked

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Link terms to ontologies.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        terms = json.load(f)

    linked = link_terms_to_ontology(terms)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(linked, f, indent=2)
