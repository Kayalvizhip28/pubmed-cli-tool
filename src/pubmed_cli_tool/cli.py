# ✅ CLI Command Handler
import click
from .core import extract_details, save_to_csv
from Bio import Entrez

@click.command()
@click.argument("query")
@click.option("--max-results", default=5, help="Number of articles to fetch")
@click.option("--file", default=None, help="CSV file path for output")
@click.option("--debug", is_flag=True, help="Enable verbose logging")
def search(query: str, max_results: int, file: str, debug: bool) -> None:
    if not query.strip():
        click.echo("❌ Error: Search query is required.")
        return

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
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

    if file and results:
        save_to_csv(file, results)
        click.echo(f"\n✅ Results saved to {file}")

if __name__ == "__main__":
    search()
