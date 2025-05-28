import json
import requests
from snakemake import input, output

with open(input[0], "r") as f:
    terms = json.load(f)

linked = {}
for term in terms:
    try:
        r = requests.get(f"https://www.ebi.ac.uk/ols/api/search?q={term}&type=class")
        data = r.json()
        if data['response']['docs']:
            linked[term] = data['response']['docs'][0].get('iri', 'Not Found')
        else:
            linked[term] = 'Not Found'
    except Exception:
        linked[term] = 'Error'

with open(output[0], "w") as f:
    json.dump(linked, f, indent=2)
