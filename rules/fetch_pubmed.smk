rule fetch_abstracts:
    output:
        "results/abstracts.json"
    params:
        email=config["email"],
        search_term=config["search_term"]
    script:
        "scripts/fetch_pubmed.py"
