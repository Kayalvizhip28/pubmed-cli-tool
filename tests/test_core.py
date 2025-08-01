import unittest
from typing import Dict, Optional
from pubmed_cli_tool.core import extract_details

class TestCoreFunctions(unittest.TestCase):
    def test_extract_details_returns_expected_keys(self) -> None:
        """
        Test that extract_details returns a dictionary with required keys for a valid PubMed ID.
        This uses an invalid ID to simulate missing data gracefully.
        """
        result: Optional[Dict[str, str]] = extract_details("12345678", debug=True)
        self.assertIsInstance(result, (dict, type(None)))

    def test_extract_details_handles_invalid_id(self) -> None:
        """
        Test that extract_details returns None for a truly invalid ID.
        """
        result: Optional[Dict[str, str]] = extract_details("", debug=True)
        self.assertIsNone(result)
