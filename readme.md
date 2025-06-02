# 🧬 Terminology Extraction Pipeline

The **Terminology Extraction Pipeline** is a modular and reproducible workflow designed for mining biomedical literature from [PubMed](https://pubmed.ncbi.nlm.nih.gov/). Powered by **Snakemake**, it automates the entire process—from abstract retrieval and term extraction to ontology linking, embedding comparison, and clustering—making it ideal for research and analysis at scale.

---

## ✨ Features

* 🔍 **PubMed Search**: Dynamically fetch abstracts using a user-defined query and email (NCBI-compliant).
* 🧠 **Terminology Extraction**: Identify biomedical terms from abstracts.
* 🔗 **Ontology Linking**: Map terms to ontologies using the [OLS API](https://www.ebi.ac.uk/ols/index).
* 🤖 **Embedding Comparison**: Compare semantic similarity using HuggingFace models like **BioBERT**, **SciBERT**, etc.
* 📊 **Clustering**: Group related terms using **DBSCAN**.
* 🛠️ **Snakemake Workflow**: Reproducible, modular, and scalable pipeline with dynamic input via command line.
* 🐳 **Docker Support**: Easily deploy the pipeline in a containerized environment.

---

## 🚀 Getting Started

### 🔧 Option 1: Local Setup (with Conda)

```bash
conda env create -f environment.yml
conda activate pubmed-env
```

Run the pipeline:

```bash
snakemake --cores 4 --config email="your@email.com" search_term="cancer biomarkers"
```

### 🐳 Option 2: Docker Setup

```bash
docker-compose up --build
```

---

## 🗂️ Project Structure

```
terminology-pipeline/
├── Snakefile
├── config.yaml
├── rules/
│   ├── fetch.smk
│   ├── extract.smk
│   ├── link.smk
│   ├── embed.smk
│   └── cluster.smk
├── scripts/
│   └── *.py
├── environment.yml
├── Dockerfile
├── docker-compose.yml
└── results/
```

---

## ⚙️ Configuration

While `email` and `search_term` are provided via the CLI, defaults and model settings live in `config.yaml`:

```yaml
output_dir: "results/"
huggingface_models:
  - dmis-lab/biobert-base-cased-v1.1
  - allenai/scibert_scivocab_uncased
```

---

## 📄 Output Files

All outputs are stored in the `results/` directory:

* `abstracts.json` – Raw PubMed abstracts
* `extracted_terms.json` – Biomedical terms extracted
* `linked_entities.json` – Ontology mappings via OLS
* `embeddings_comparison.csv` – Embedding similarity scores
* `clusters.json` – Grouped terms from DBSCAN

---

## 📦 Dependencies

Dependencies are handled via Conda. Defined in `environment.yml`:

```yaml
- biopython
- pandas
- tqdm
- scikit-learn
- requests
- snakemake
- sentence-transformers (via pip)
```

---

## 🙋 Contributing

Open to pull requests and discussions. Feel free to contribute enhancements or report issues!

---

## 📄 License

MIT License. See [`LICENSE`](LICENSE) for more information.

---


