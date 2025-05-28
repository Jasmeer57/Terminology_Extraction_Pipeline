import json
from Bio import Entrez
from snakemake import params, output

Entrez.email = params.email
search_term = params.search_term

handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100)
record = Entrez.read(handle)
ids = record["IdList"]

abstracts = []
for pubmed_id in ids:
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="abstract", retmode="text")
    abstract_text = handle.read()
    abstracts.append({"id": pubmed_id, "abstract": abstract_text})

with open(output[0], "w") as f:
    json.dump(abstracts, f, indent=2)
