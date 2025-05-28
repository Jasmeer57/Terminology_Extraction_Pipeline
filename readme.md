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

---

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/pubmed-pipeline-tool.git
cd pubmed-pipeline-tool
```

### 2. Install Python Dependencies

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

---

## 🐳 Docker Deployment

### Using Docker Compose (Recommended)

```bash
docker-compose up --build
```

### Using Docker Only

```bash
docker build -t pubmed-pipeline .
docker run -p 8501:8501 pubmed-pipeline
```

The GUI will be accessible at `http://localhost:8501`.

---

## 📂 Output

The tool generates outputs in the specified directory, including:

* `abstracts.json` — Raw PubMed abstracts.
* `extracted_terms.json` — Extracted biomedical terms.
* `linked_entities.json` — Terms linked to ontology concepts via OLS API.
* `embeddings_comparison.csv` — Similarity and comparison scores across multiple HuggingFace models.
* `clusters.json` — DBSCAN-generated clusters of related terms.

---

## 🧪 Example Use Case

1. Launch the GUI (locally or via Docker).
2. Enter your PubMed search query (e.g., *"Alzheimer's disease biomarkers"*).
3. Input your email for NCBI compliance.
4. Click **Run** and wait for the processing to complete.
5. Explore and utilize the structured output data for research, machine learning, or visualization tasks.

---

## 🛠️ Development

To run the tool locally without Docker:

```bash
pip install -r requirements.txt
streamlit run app.py
```

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
