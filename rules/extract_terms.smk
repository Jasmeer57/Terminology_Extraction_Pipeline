rule extract_terms:
    input:
        "results/abstracts.json"
    output:
        "results/extracted_terms.json"
    script:
        "scripts/extract_terms.py"
