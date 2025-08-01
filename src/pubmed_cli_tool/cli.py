# src/pubmed_cli_tool/cli.py

import click
from typing import Optional
from Bio import Entrez
from pubmed_cli_tool.core import extract_details, save_to_csv


@click.command()
@click.argument("query")
@click.option(
    "--max-results",
    default=5,
    type=int,
    help="Maximum number of articles to fetch from PubMed (default is 5)."
)
@click.option(
    "--file",
    default=None,
    type=str,
    help="CSV file path to save the filtered output."
)
@click.option(
    "--debug",
    is_flag=True,
    help="Enable debug logging to display extra processing information."
)
def search(query: str, max_results: int, file: Optional[str], debug: bool) -> None:
    """
    CLI command to search PubMed for articles matching the query,
    extract industry-affiliated article metadata, and optionally save to CSV.

    Args:
        query (str): Search query for PubMed.
        max_results (int): Maximum number of articles to fetch.
        file (Optional[str]): Path to save output CSV file.
        debug (bool): Enable debug logging if True.
    """
    if debug:
        click.echo(f"🔍 Searching PubMed for: {query}, max results: {max_results}")

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])
    except Exception as e:
        click.echo(f"❌ PubMed search failed: {e}")
        return

    if not ids:
        click.echo("⚠️ No results found.")
        return

    results = []
    for pubmed_id in ids:
        details = extract_details(pubmed_id, debug)
        if details:
            results.append(details)
            click.echo("\n--- Article ---")
            for key, value in details.items():
                click.echo(f"{key}: {value}")
        elif debug:
            click.echo(f"⚠️ No industry affiliation found for PubMed ID: {pubmed_id}")

    if file and results:
        try:
            save_to_csv(file, results)
            click.echo(f"\n✅ Results saved to {file}")
        except Exception as e:
            click.echo(f"❌ Failed to save CSV file: {e}")
    elif not results:
        click.echo("⚠️ No industry-related articles found.")
