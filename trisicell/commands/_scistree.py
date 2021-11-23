import os

import click

import trisicell as tsc


@click.command(short_help="Run ScisTree.")
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
    "--n_threads",
    "-p",
    default=1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def scistree(genotype_file, alpha, beta, n_threads):
    """ScisTree.

    Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach
    :cite:`ScisTree`.

    trisicell scistree input.SC 0.0001 0.1 -p 1
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.scistree.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.scistree(df_in, alpha=alpha, beta=beta, n_threads=n_threads)
    tsc.io.write(df_out, f"{outfile}.scistree.CFMatrix")

    return None
