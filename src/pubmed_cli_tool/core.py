# ✅ Core Logic File
import os, csv, logging
from typing import Optional, Dict, List
from dotenv import load_dotenv
from Bio import Entrez

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "PubMedCLI"
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def is_industry_affiliation(affiliation: str) -> bool:
    keywords = ["Inc", "Ltd", "Pharma", "Biotech", "Therapeutics", "Diagnostics", "Corp", "LLC"]
    return any(kw.lower() in affiliation.lower() for kw in keywords)

def extract_details(pubmed_id: str, debug: bool = False) -> Optional[Dict[str, str]]:
    # ✅ Copied from your working logic
    # (Make sure it's complete as per your previous version)
    ...

def save_to_csv(filename: str, articles: List[Dict[str, str]]) -> None:
    # ✅ Copied logic here too
    ...
