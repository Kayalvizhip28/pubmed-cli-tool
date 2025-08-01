import unittest
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from pubmed_cli_tool.cli import search

class TestPubMedCLI(unittest.TestCase):
    def setUp(self) -> None:
        self.runner = CliRunner()

    @patch("pubmed_cli_tool.cli.extract_details")
    @patch("Bio.Entrez.read")
    @patch("Bio.Entrez.esearch")
    def test_valid_query_runs(self, mock_esearch, mock_read, mock_extract):
        """
        Simulate a valid PubMed search that returns an industry-related article.
        """
        mock_esearch.return_value = MagicMock()
        mock_read.return_value = {"IdList": ["12345678"]}
        mock_extract.return_value = {
            "PubMed ID": "12345678",
            "Title": "Sample Article",
            "Publication Date": "2024-06",
            "Non-academic Author(s)": "John Doe",
            "Company Affiliations": "BioTech Inc",
            "Corresponding Email": "sample@example.com"
        }

        result = self.runner.invoke(
            search, ["cancer", "--max-results", "1", "--debug"]
        )

        self.assertEqual(result.exit_code, 0)
        self.assertIn("Title: Sample Article", result.output)
