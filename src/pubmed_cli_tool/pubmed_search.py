import os
import csv
import logging
import click
from typing import Optional, Dict, List
from dotenv import load_dotenv
from Bio import Entrez

# Load environment variables
load_dotenv()

# Set up Entrez identity
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "PubMedCLI"

# Configure logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def is_industry_affiliation(affiliation: str) -> bool:
    """
    Check if the affiliation string suggests an industry-based institution.

    Args:
        affiliation (str): Author's affiliation description

    Returns:
        bool: True if affiliation contains industry keywords
    """
    keywords = ["Inc", "Ltd", "Pharma", "Biotech", "Therapeutics", "Diagnostics", "Corp", "LLC"]
    return any(kw.lower() in affiliation.lower() for kw in keywords)

def extract_details(pubmed_id: str, debug: bool = False) -> Optional[Dict[str, str]]:
    """
    Fetch and extract article metadata for a given PubMed ID.

    Args:
        pubmed_id (str): Unique identifier of the PubMed article
        debug (bool): Enables debug logging

    Returns:
        Optional[Dict[str, str]]: Article details if industry-affiliated author is found; None otherwise
    """
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        records = Entrez.read(handle)
        article = records["PubmedArticle"][0]
    except Exception as e:
        if debug:
            logger.error(f"Error fetching details for PubMed ID {pubmed_id}: {e}")
        return None

    # Extract title and publication date
    try:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
        authors = article["MedlineCitation"]["Article"].get("AuthorList", [])
    except KeyError as e:
        if debug:
            logger.warning(f"Missing expected fields in article {pubmed_id}: {e}")
        return None

    non_academic: List[str] = []
    companies: List[str] = []
    email: Optional[str] = None
    found_industry = False

    for author in authors:
        affiliation = author.get("AffiliationInfo", [{}])[0].get("Affiliation", "")
        if is_industry_affiliation(affiliation):
            found_industry = True
            full_name = f"{author.get('ForeName', '')} {author.get('LastName', '')}".strip()
            non_academic.append(full_name)
            companies.append(affiliation)

            # Extract first encountered email (if any)
            if "@" in affiliation and not email:
                for word in affiliation.split():
                    if "@" in word:
                        email = word.strip("<>.,;")

    if not found_industry:
        return None

    pubdate = f"{pub_date.get('Year', '')}-{pub_date.get('Month', '')}"

    return {
        "PubMed ID": pubmed_id,
        "Title": title,
        "Publication Date": pubdate,
        "Non-academic Author(s)": "; ".join(non_academic),
        "Company Affiliations": "; ".join(companies),
        "Corresponding Email": email or "N/A"
    }

def save_to_csv(filename: str, articles: List[Dict[str, str]]) -> None:
    """
    Write article results to a CSV file.

    Args:
        filename (str): Destination file path
        articles (List[Dict[str, str]]): List of article metadata
    """
    try:
        with open(filename, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=articles[0].keys())
            writer.writeheader()
            writer.writerows(articles)
    except Exception as e:
        logger.error(f"Failed to save CSV: {e}")

@click.command()
@click.argument("query")
@click.option("--max-results", default=5, help="Number of articles to return")
@click.option("--file", default=None, help="Output CSV file path")
@click.option("--debug", is_flag=True, help="Enable debug logging")
def search(query: str, max_results: int, file: Optional[str], debug: bool) -> None:
    """
    CLI entry point to search PubMed and filter articles with industry authors.

    Args:
        query (str): Search query string
        max_results (int): Number of articles to retrieve
        file (Optional[str]): CSV file to export results
        debug (bool): Flag to enable verbose logging
    """
    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug mode activated")

    # Basic validation for empty search query
    if not query.strip():
        click.echo("❌ Error: Search query cannot be empty.")
        return

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
    except Exception as e:
        click.echo(f"❌ Failed to retrieve articles: {e}")
        return

    if not ids:
        click.echo("⚠️ No articles found for the given query.")
        return

    results: List[Dict[str, str]] = []
    for pubmed_id in ids:
        details = extract_details(pubmed_id, debug)
        if details:
            results.append(details)
            click.echo("\n--- Article ---")
            for key, value in details.items():
                click.echo(f"{key}: {value}")

    if file and results:
        save_to_csv(file, results)
        click.echo(f"\n✅ Results saved to '{file}'")

if __name__ == "__main__":
    search()

