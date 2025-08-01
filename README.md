# 🧬PubMed CLI Tool — `get-papers-list`

A Python-based command-line interface (CLI) to search PubMed articles and filter authors with **non-academic (industry)** affiliations such as Biotech, Pharma, Diagnostics, etc.

This tool fetches relevant articles using the **NCBI Entrez API**, parses metadata like author names, affiliations, publication date, and email, and saves the results to a CSV file.

---

## Project Structure

pubmed_cli_tool/
├── src/
│ └── pubmed_cli_tool/
│ ├── cli.py # CLI command handler using Click
│ └── core.py # Business logic for data fetching, filtering, and CSV saving
├── tests/
│ ├── test_core.py # Unit tests for core logic
│ └── test_cli.py # Unit tests for CLI interactions using CliRunner
├── .env # Contains sensitive config (e.g., NCBI_EMAIL)
├── pyproject.toml # Poetry project config and CLI entry definition
└── README.md # Project documentation (you are here!)


---

## Installation (via Poetry)
Poetry will automatically create a virtual environment and install all required packages:

poetry install

## Clone the Repository
git clone https://github.com/your-username/pubmed_cli_tool.git

cd pubmed_cli_tool

Replace your-username with your actual GitHub username.

## Environment Setup
Create a .env file in the project root with your registered email to use NCBI’s Entrez API.

NCBI_EMAIL=your.email@example.com

This email is required by NCBI to monitor usage and prevent abuse.

## Usage
Use the command-line tool to search PubMed and extract article metadata.

Example

poetry run get-papers-list "cancer immunotherapy" --max-results 3 --file results.csv --debug

This command will:

Fetch top 3 articles matching "cancer immunotherapy"

Filter authors with industry affiliations

Save results to results.csv

Print debug information if enabled

## Optional Flags
Option	Description

--file	Save the results to a CSV file

--max-results	Number of articles to fetch (default: 5)

--debug	 Enable verbose logging and error tracing

## Exposing CLI as get-papers-list
Already configured in pyproject.toml:

[tool.poetry.scripts]

get-papers-list = "pubmed_cli_tool.cli:search"

Now you can run the tool like this:

poetry run get-papers-list "breast cancer" --max-results 2

## Libraries Used
- [Click](https://click.palletsprojects.com/) – elegant CLI creation
- [Biopython](https://biopython.org/) – access to NCBI Entrez API
- [python-dotenv](https://github.com/theskumar/python-dotenv) – manage environment variables
- [Poetry](https://python-poetry.org/) – dependency management & packaging
- [unittest](https://docs.python.org/3/library/unittest.html) – built-in Python testing framework

## Running Tests
poetry run pytest

This runs both:

tests/test_core.py – tests for business logic

tests/test_cli.py – tests for CLI behavior using Click’s CliRunner

## CSV Output Format
Example CSV output (results.csv):

PubMed ID	  Title	   Publication Date	   Non-academic Author(s)	  Company Affiliations	  Corresponding Email

40739319	  In-vitro assessment...	  2025-Jul	   Rupali Ghosh; Noor Fatima	J  amia Hamdard Biotech Lab	   swajid@jamiahamdard.ac.in

Each row represents a filtered article with at least one industry-affiliated author.




