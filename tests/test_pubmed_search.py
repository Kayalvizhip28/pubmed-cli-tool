import unittest
from unittest.mock import patch, MagicMock
from src.pubmed_search import extract_details
from Bio import Entrez

class TestPubMedSearch(unittest.TestCase):

    @patch("src.pubmed_search.Entrez.efetch")
    @patch("src.pubmed_search.Entrez.read")
    def test_extract_details_no_industry_affiliation(self, mock_read, mock_efetch):
        # Simulated XML response from Entrez.read with only academic affiliation
        mock_article = {
            "PubmedArticle": [{
                "MedlineCitation": {
                    "Article": {
                        "ArticleTitle": "Academic Research on Cancer",
                        "Journal": {
                            "JournalIssue": {
                                "PubDate": {
                                    "Year": "2023",
                                    "Month": "07"
                                }
                            }
                        },
                        "AuthorList": [{
                            "ForeName": "Alice",
                            "LastName": "Smith",
                            "AffiliationInfo": [{
                                "Affiliation": "Department of Oncology, University of Somewhere"
                            }]
                        }]
                    }
                }
            }]
        }

        mock_efetch.return_value = MagicMock()
        mock_read.return_value = mock_article

        result = extract_details("12345678")
        self.assertIsNone(result, "Should return None for academic-only affiliations")

if __name__ == "__main__":
    unittest.main()
