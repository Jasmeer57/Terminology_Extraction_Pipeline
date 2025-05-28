import json
from snakemake import input, output

with open(input[0], "r") as f:
    abstracts = json.load(f)

# Naive term extraction (for demo purposes)
terms = set()
for entry in abstracts:
    for word in entry["abstract"].split():
        if word[0].isupper() and len(word) > 3:
            terms.add(word.strip(",.;"))

with open(output[0], "w") as f:
    json.dump(list(terms), f, indent=2)
