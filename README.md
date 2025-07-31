# 🔎 PubMed CLI Tool

A Python-based command-line interface to search **PubMed** articles and extract results involving **non-academic authors**—such as those affiliated with biotech or pharmaceutical companies. Ideal for filtering publications with industry relevance and exporting them for further analysis.

---

## 🗂️ Project Structure

```plaintext
pubmed_cli_tool/
├── src/
│   └── pubmed_search.py        # CLI logic and article filtering
├── tests/
│   └── test_pubmed_search.py   # Unit tests with mock Entrez data
├── .env                        # Stores NCBI_EMAIL environment variable
├── README.md                   # Project documentation
├── pyproject.toml              # Poetry dependency configuration
