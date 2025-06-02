rule fetch_pubmed:
    output:
        "results/abstracts.json"
    params:
        email=config["email"],
        term=config["search_term"]
    shell:
        """
        python scripts/fetch_pubmed.py --email {params.email} --term "{params.term}" --out {output}
        """
