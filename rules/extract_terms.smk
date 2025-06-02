rule extract_terms:
    input:
        "results/abstracts.json"
    output:
        "results/extracted_terms.json"
    shell:
        """
        python scripts/extract_terms.py --input {input} --output {output}
        """
