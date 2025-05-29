# 🧬 Terminology Extraction Pipeline Tool (PubMed Pipeline Tool)

The **Terminology Extraction Pipeline Tool** is a comprehensive, user-friendly application designed to streamline biomedical research workflows by automating the process of querying [PubMed](https://pubmed.ncbi.nlm.nih.gov/), extracting abstracts, performing terminology extraction, entity linking to biomedical ontologies, and analyzing term embeddings using state-of-the-art models — all accessible via an intuitive graphical user interface (GUI). The tool also supports Docker for easy deployment across environments.

---

## ✨ Features

* 🔍 **PubMed Search**
  Query PubMed with user-provided search strings and retrieve relevant abstracts.

* 📝 **Abstract Extraction**
  Automatically fetch abstracts from PubMed articles for downstream analysis.

* 🧠 **Term Extraction**
  Extract biomedical terms from abstracts to identify key concepts.

* 🔗 **Entity Linking**
  Link extracted terms to ontologies using the [Ontology Lookup Service (OLS) API](https://www.ebi.ac.uk/ols/index) to provide semantic context.

* 🤖 **Model Comparison**
  Generate and compare embeddings of extracted terms using multiple HuggingFace transformer models such as `BioBERT`, `SciBERT`, and others.

* 📊 **Clustering**
  Group terms via the DBSCAN clustering algorithm to identify meaningful clusters and relationships.

* 🖥️ **Intuitive GUI**
  A simple graphical interface allows easy input of search queries and user credentials, streamlining interaction.

* 🐳 **Docker Support**
  Includes pre-configured Dockerfile and `docker-compose.yml` for hassle-free containerized deployment.

* ⚙️ **Customizable Configuration**
  Configure parameters like search term, output directory, and model choices through a `config.json` file.
  
* 🛠️ **Snakemake Workflow**
  Reproducible, modular, and scalable pipeline with dynamic input via command line.

---

## 🚀 Getting Started

### Download Guide
#### 1. Clone the Repository

```bash
git clone https://github.com/your-username/pubmed-pipeline-tool.git
cd pubmed-pipeline-tool
```

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


### Option 3. Install Python Dependencies

```bash
pip install biopython pandas tqdm transformers sentence-transformers scikit-learn requests tkinter
```

Alternatively, all dependencies are listed in `requirements.txt` and will be installed automatically when using Docker.

### 3. Configure

Edit the `config.json` file to set your email, search term, output directory, and preferred HuggingFace models:

```json
{
  "email": "your-email@example.com",
  "search_term": "cancer genomics",
  "output_dir": "./results",
  "huggingface_models": [
    "dmis-lab/biobert-base-cased-v1.1",
    "allenai/scibert_scivocab_uncased"
  ]
}
```
The GUI will be accessible at `http://localhost:8501`.


## 🧪 Example Use Case

1. Launch the GUI (locally or via Docker).
2. Enter your PubMed search query (e.g., *"Alzheimer's disease biomarkers"*).
3. Input your email for NCBI compliance.
4. Click **Run** and wait for the processing to complete.
5. Explore and utilize the structured output data for research, machine learning, or visualization tasks.


---

## 🙋‍♀️ Contributing

Contributions, issues, and feature requests are welcome! Feel free to open an issue or submit a pull request.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

## 📬 Contact

For questions or collaboration, please contact:

**Jasmeer Singh Kalra**
📧 [jasmeer57@gmail.com](mailto:jasmeer57@gmail.com)
🔗 [LinkedIn](https://www.linkedin.com/in/jasmeer-singh/)
🐙 [GitHub](https://github.com/Jasmeer57)

---

## ⭐ Acknowledgments

* HuggingFace for state-of-the-art biomedical transformer models
* NCBI and PubMed for access to invaluable biomedical research data
* EBI Ontology Lookup Service (OLS) for ontology linking capabilities

---
