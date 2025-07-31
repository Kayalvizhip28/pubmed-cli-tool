import os
import csv
import click
import logging
from dotenv import load_dotenv
from Bio import Entrez

# Load environment variables
load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "PubMedCLI"

# Setup logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def is_industry_affiliation(affiliation):
    keywords = ["Inc", "Ltd", "Pharma", "Biotech", "Therapeutics", "Diagnostics", "Corp", "LLC"]
    return any(kw.lower() in affiliation.lower() for kw in keywords)

def extract_details(pubmed_id, debug=False):
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        records = Entrez.read(handle)
        article = records["PubmedArticle"][0]
    except Exception as e:
        if debug:
            logger.warning(f"Failed to fetch details for {pubmed_id}: {e}")
        return None

    title = article["MedlineCitation"]["Article"]["ArticleTitle"]
    pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
    authors = article["MedlineCitation"]["Article"].get("AuthorList", [])

    non_academic = []
    companies = []
    email = None
    found_industry = False

    for author in authors:
        affiliation = author.get("AffiliationInfo", [{}])[0].get("Affiliation", "")
        if is_industry_affiliation(affiliation):
            found_industry = True
            full_name = f"{author.get('ForeName', '')} {author.get('LastName', '')}".strip()
            non_academic.append(full_name)
            companies.append(affiliation)
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

def save_to_csv(filename, articles):
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=articles[0].keys())
        writer.writeheader()
        writer.writerows(articles)

@click.command()
@click.argument("query")
@click.option("--max-results", default=5, help="Number of articles to return")
@click.option("--file", default=None, help="Output CSV file path")
@click.option("--debug", is_flag=True, help="Enable debug logging")
def search(query, max_results, file, debug):
    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug mode activated")

    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]

    results = []
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
