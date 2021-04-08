import os

import click

import trisicell as tsc


@click.command(short_help="Run HUNTRESS.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "alpha",
    required=True,
    type=float,
)
@click.argument(
    "beta",
    required=True,
    type=float,
)
@click.option(
    "--method",
    "-m",
    default="dna",
    type=click.Choice(["dna", "fn", "rna"]),
    show_default=True,
    help="Method of the HUNTRESS",
)
@click.option(
    "--n_threads",
    "-p",
    default=1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def huntress(genotype_file, alpha, beta, method, n_threads):
    """HUNTRESS

    trisicell huntress input.SC 0.0001 0.1 -m dna -p 8
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.huntress.log"

    df_out = tsc.tl.huntress(
        genotype_file, alpha=alpha, beta=beta, kind=method, n_threads=n_threads
    )
    tsc.io.write(df_out, f"{outfile}.huntress.CFMatrix")

    return None
