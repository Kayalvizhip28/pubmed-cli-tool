# 🔎 PubMed CLI Tool

A Python-based command-line interface to search **PubMed** articles and extract metadata for papers involving **non-academic authors**—such as those affiliated with **biotech, pharma, diagnostics, or industry-related organizations**. This tool is ideal for research analysts, medical scientists, and industry professionals to extract filtered research papers and export results into CSV format for further analysis.

---

## 🗂️ Project Structure

```plaintext
pubmed_cli_tool/
├── src/
│   └── pubmed_cli_tool/
│       └── pubmed_search.py       # Core CLI logic: fetches, filters, and exports results
├── tests/
│   └── test_pubmed_search.py      # Unit tests with mocked Entrez API responses
├── .env                           # Stores NCBI_EMAIL environment variable (ignored by Git)
├── README.md                      # Project documentation (you are here)
├── pyproject.toml                 # Poetry dependency and script configuration
└── .gitignore                     # Files/directories excluded from version control

📦 Installation & Setup
✅ Requires Python 3.10+ and Poetry for dependency management.

1. Clone the repository:

bash
Copy
Edit
git clone https://github.com/your-username/pubmed_cli_tool.git
cd pubmed_cli_tool

2. Install dependencies using Poetry:

bash
Copy
Edit
poetry install

3. Set up environment variable for NCBI Entrez access:

Create a .env file in the project root:

ini
Copy
Edit
NCBI_EMAIL=your_email@example.com
This email is used by NCBI to identify your API usage.

🚀 Usage
Run the CLI tool using Poetry and provide a search query (in quotes):

bash
Copy
Edit
poetry run get-papers-list "cancer immunotherapy"

🔧 Optional Flags
Flag	Description	Example
--max-results	Max number of PubMed articles to fetch (default: 5)	--max-results 10
--file	Save results to a CSV file	--file cancer_results.csv
--debug	Enable verbose logging for troubleshooting	--debug

✅ Example
bash
Copy
Edit
poetry run get-papers-list "cancer immunotherapy" --max-results 10 --file results.csv --debug
This command:

Searches PubMed for up to 10 papers matching the query

Filters for non-academic/industry-affiliated authors

Saves the extracted metadata into results.csv

Displays logs in debug mode

🛠️ Tools & Libraries Used
🧬 Biopython – Used for interacting with NCBI Entrez APIs

⚙️ Click – For creating clean and flexible CLI commands

🧪 python-dotenv – For securely managing environment variables

📦 Poetry – Dependency and package management

🧪 Built-in unittest.mock – For mocking API calls during testing

🧪 Running Tests
Run unit tests to validate CLI behavior:

bash
Copy
Edit
poetry run pytest
Tests are located in the /tests folder and mock the Entrez API using test cases like:

Valid and invalid PubMed IDs

Articles with no non-academic affiliations

Email detection edge cases

⚙️ Setting Up Executable CLI Command
Your CLI command is exposed as get-papers-list via Poetry.

Inside pyproject.toml:

toml
Copy
Edit
[tool.poetry.scripts]
get-papers-list = "pubmed_cli_tool.pubmed_search:search"
✅ This lets you run the tool like:

bash
Copy
Edit
poetry run get-papers-list "your query here"
No need to directly call the Python script!

🧹 .gitignore Configuration
To avoid pushing sensitive and temporary files to Git:

markdown
Copy
Edit
# .gitignore
.env
__pycache__/
*.py[cod]
*.csv
Make sure to commit your .gitignore after editing:

bash
Copy
Edit
git add .gitignore
git commit -m "Update .gitignore to exclude secrets and build artifacts"
📤 Publishing to GitHub
Initialize and push your project:

bash
Copy
Edit
git init
git remote add origin https://github.com/your-username/pubmed_cli_tool.git
git add .
git commit -m "Initial commit"
git push -u origin main
✅ Confirm your code is published at:

arduino
Copy
Edit
https://github.com/your-username/pubmed_cli_tool
📎 Output Example (CSV)
A sample row from the generated CSV:

csv
Copy
Edit
PubMed ID,Title,Publication Date,Non-academic Author(s),Company Affiliations,Corresponding Email
40738546,"Nanoparticles for mRNA-based cancer immunotherapy.","2025-Jul","Rakesh Pahwa; Gulshan Sharma","Institute