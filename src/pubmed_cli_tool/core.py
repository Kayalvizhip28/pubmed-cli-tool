import os
import csv
import logging
from typing import Optional, Dict, List
from dotenv import load_dotenv
from Bio import Entrez

# Load environment variables from .env
load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "PubMedCLI"

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def is_industry_affiliation(affiliation: str) -> bool:
    """
    Check if the affiliation contains industry-related keywords.
    """
    keywords = ["Inc", "Ltd", "Pharma", "Biotech", "Therapeutics", "Diagnostics", "Corp", "LLC"]
    return any(kw.lower() in affiliation.lower() for kw in keywords)

def extract_details(pubmed_id: str, debug: bool = False) -> Optional[Dict[str, str]]:
    """
    Fetch and extract industry-relevant metadata for a given PubMed ID.
    Returns None if no industry affiliation is found.
    """
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        article = records["PubmedArticle"][0]
    except Exception as e:
        if debug:
            logger.error(f"❌ Error fetching details for PubMed ID {pubmed_id}: {e}")
        return None

    try:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
        authors = article["MedlineCitation"]["Article"].get("AuthorList", [])
    except KeyError as e:
        if debug:
            logger.warning(f"⚠️ Missing expected fields in article {pubmed_id}: {e}")
        return None

    non_academic: List[str] = []
    companies: List[str] = []
    email: Optional[str] = None
    found_industry = False

    for author in authors:
        aff_info = author.get("AffiliationInfo", [])
        if not aff_info:
            continue

        affiliation = aff_info[0].get("Affiliation", "")
        if is_industry_affiliation(affiliation):
            found_industry = True
            full_name = f"{author.get('ForeName', '')} {author.get('LastName', '')}".strip()
            non_academic.append(full_name)
            companies.append(affiliation)

            # Extract email if not already found
            if "@" in affiliation and not email:
                for word in affiliation.split():
                    if "@" in word:
                        email = word.strip("<>.,;()[]")

    if not found_industry:
        if debug:
            logger.info(f"⚠️ No industry affiliation found for PubMed ID: {pubmed_id}")
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
    Save a list of article metadata dictionaries to a CSV file.
    """
    try:
        with open(filename, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=articles[0].keys())
            writer.writeheader()
            writer.writerows(articles)
    except Exception as e:
        logger.error(f"❌ Failed to save CSV to {filename}: {e}")
