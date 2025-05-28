rule entity_linking:
    input:
        "results/extracted_terms.json"
    output:
        "results/linked_entities.json"
    script:
        "scripts/link_entities.py"
