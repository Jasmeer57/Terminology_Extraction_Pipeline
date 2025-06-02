rule entity_linking:
    input:
        "results/extracted_terms.json"
    output:
        "results/linked_entities.json"
    shell:
        """
        python scripts/entity_linking.py --input {input} --output {output}
        """
