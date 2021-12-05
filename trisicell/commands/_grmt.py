import os

import click

import trisicell as tsc


@click.command(short_help="Run GRMT.")
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
    "--n_iters",
    "-l",
    default=30,
    type=int,
    show_default=True,
    help="Number of iterations.",
)
@click.option(
    "--n_threads",
    "-p",
    default=1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def grmt(genotype_file, alpha, beta, n_iters, n_threads):
    """GRMT.

    Generative Reconstruction of Mutation Tree From Scratch Using Single-Cell
    Sequencing Data :cite:`GRMT`.

    trisicell grmt input.SC 0.0001 0.1 -l 100 -p 2
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.grmt.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.grmt(
        df_in, alpha=alpha, beta=beta, n_iters=n_iters, n_threads=n_threads
    )
    tsc.io.write(df_out, f"{outfile}.grmt.CFMatrix")

    return None
