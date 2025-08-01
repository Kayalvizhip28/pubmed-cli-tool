import click
from typing import Optional
from Bio import Entrez
from pubmed_cli_tool.core import extract_details, save_to_csv

@click.command()
@click.argument("query")
@click.option("--max-results", default=5, help="Maximum number of articles to fetch.")
@click.option("--file", default=None, help="CSV file path for output.")
@click.option("--debug", is_flag=True, help="Enable verbose logging.")
def search(query: str, max_results: int, file: Optional[str], debug: bool) -> None:
    """
    CLI command to search PubMed articles and optionally save filtered results to CSV.
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
        save_to_csv(file, results)
        click.echo(f"\n✅ Results saved to {file}")
    elif not results:
        click.echo("⚠️ No industry-related articles found.")
