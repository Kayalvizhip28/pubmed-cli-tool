import unittest
from typing import Optional, Dict
from pubmed_cli_tool.core import extract_details


class TestCoreFunctions(unittest.TestCase):
    """
    Unit tests for core logic in pubmed_cli_tool.core module.
    Verifies data structure, edge case handling, and simulated behavior of extract_details().
    """

    def test_extract_details_returns_expected_keys(self) -> None:
        """
        Test that extract_details returns a dictionary with required keys for a valid PubMed ID.
        This uses the dummy implementation that returns a sample article.
        """
        result: Optional[Dict[str, str]] = extract_details("12345678")
        self.assertIsInstance(result, dict)
        self.assertIn("PubMed ID", result)
        self.assertIn("Title", result)
        self.assertIn("Corresponding Email", result)

    def test_extract_details_handles_invalid_id(self) -> None:
        """
        Test that extract_details still returns a valid structure for an empty or invalid ID.
        This simulates graceful fallback behavior for malformed input.
        """
        result: Optional[Dict[str, str]] = extract_details("")
        self.assertIsNotNone(result)
        self.assertEqual(result.get("PubMed ID", ""), "")

if __name__ == "__main__":
    unittest.main()
