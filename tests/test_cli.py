import unittest
from click.testing import CliRunner
from pubmed_cli_tool.cli import search
from unittest.mock import patch, MagicMock


class TestPubMedCLI(unittest.TestCase):
    """
    Unit tests for the PubMed CLI search tool using Click's CliRunner.
    Covers CLI argument validation, search result handling, and mocking API interactions.
    """

    def setUp(self) -> None:
        """Initialize CLI test runner before each test."""
        self.runner: CliRunner = CliRunner()

    def test_missing_query_argument(self) -> None:
        """
        Test that the CLI returns an error when no query argument is provided.
        """
        result = self.runner.invoke(search, [])
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Error: Missing argument 'QUERY'", result.output)

    @patch("pubmed_cli_tool.core.extract_details")
    @patch("Bio.Entrez.read")
    @patch("Bio.Entrez.esearch")
    def test_invalid_query_returns_no_results(
        self,
        mock_esearch: MagicMock,
        mock_read: MagicMock,
        mock_extract: MagicMock,
    ) -> None:
        """
        Test a valid command with an invalid search term that yields no PubMed results.
        """
        mock_esearch.return_value = MagicMock()
        mock_read.return_value = {"IdList": []}
        mock_extract.return_value = None

        result = self.runner.invoke(search, ["invalidquery", "--max-results", "1"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("⚠️ No results found.", result.output)

    @patch("pubmed_cli_tool.core.extract_details")
    @patch("Bio.Entrez.read")
    @patch("Bio.Entrez.esearch")
    def test_valid_query_runs(
        self,
        mock_esearch: MagicMock,
        mock_read: MagicMock,
        mock_extract: MagicMock,
    ) -> None:
        """
        Test a valid search query and successful data fetch and display.
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

        result = self.runner.invoke(search, ["cancer", "--max-results", "1"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Title: Sample Article", result.output)
        self.assertIn("Corresponding Email: sample@example.com", result.output)
        self.assertIn("Company Affiliations: BioTech Inc", result.output)

if __name__ == "__main__":
    unittest.main()
