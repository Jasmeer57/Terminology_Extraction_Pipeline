import argparse
import json
from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch
from tqdm import tqdm

# Load BioBERT model
model_name = "dmis-lab/biobert-base-cased-v1.1"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name)
model.eval()

def score_term(term):
    sentence = f"The role of {term} in gut microbiome regulation has been studied."
    inputs = tokenizer(sentence, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits
        term_ids = tokenizer(term, add_special_tokens=False).input_ids
        if not term_ids:
            return 0.0  # fallback if term not tokenized
        score = sum(logits[0, i, term_id].item() for i, term_id in enumerate(term_ids)) / len(term_ids)
    return score

def extract_terms(abstracts):
    term_scores = {}
    for entry in tqdm(abstracts, desc="Extracting terms"):
        text = entry.get("Abstract", "")
        tokens = tokenizer.tokenize(text)
        # Remove special tokens
        tokens = [t for t in tokens if not t.startswith("##") and t.isalpha() and len(t) > 2]
        unique_tokens = set(tokens)
        for token in unique_tokens:
            if len(token) < 3 or len(token) > 50:
                continue
            try:
                score = score_term(token)
                if score > 5.0:  # arbitrary threshold for relevance
                    term_scores[(token.lower())] = {
                        "term": token,
                        "score": score
                    }
            except Exception as e:
                print(f"Error scoring {token}: {e}")

    # Convert to list of dicts
    terms = list(term_scores.values())
    return terms

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract terms using BioBERT scoring")
    parser.add_argument("--input", required=True, help="Input JSON file with abstracts")
    parser.add_argument("--output", required=True, help="Output JSON file with extracted terms")
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        abstracts = json.load(f)

    terms = extract_terms(abstracts)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(terms, f, indent=2, ensure_ascii=False)

    print(f"Extracted {len(terms)} unique terms and saved to {args.output}")
