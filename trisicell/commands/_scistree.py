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
def scistree(genotype_file, alpha, beta):
    """Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach :cite:`ScisTree`.

    trisicell scistree input.SC 0.0001 0.1
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.scistree.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.scistree(df_in, alpha=alpha, beta=beta)
    tsc.io.write(df_out, f"{outfile}.scistree.CFMatrix")

    return None
