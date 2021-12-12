import os

import click

import trisicell as tsc


@click.command(short_help="Run gpps.")
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
    "--k_dollo",
    "-k",
    default=0,
    type=int,
    show_default=True,
    help="k-Dollo.",
)
@click.option(
    "--max_del",
    "-d",
    default=-1,
    type=int,
    show_default=True,
    help="Maximum number of deletion allowed.",
)
@click.option(
    "--neighbor_size",
    "-s",
    default=30,
    type=int,
    show_default=True,
    help="Hill climbing neighborhood size.",
)
@click.option(
    "--n_iters",
    "-l",
    default=100,
    type=int,
    show_default=True,
    help="Hill climbing maximum iterations.",
)
@click.option(
    "--time_limit",
    "-t",
    default=86400,
    type=int,
    show_default=True,
    help="Time limit of the program (in second).",
)
@click.option(
    "--n_threads",
    "-p",
    default=1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def gpps(
    genotype_file,
    alpha,
    beta,
    k_dollo,
    max_del,
    neighbor_size,
    n_iters,
    time_limit,
    n_threads,
):
    """gpps.

    an ILP-based approach for inferring cancer progression with mutation losses from
    single cell data :cite:`gpps`.

    trisicell gpps input.SC 0.0001 0.1 -k 0 -s 30 -l 100 -t 86400 -p 1
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.gpps.log"

    df_in = tsc.io.read(genotype_file)

    if alpha == 0:
        alpha = 0.000000000001
    df_out = tsc.tl.gpps(
        df_in,
        alpha=alpha,
        beta=beta,
        k_dollo=k_dollo,
        max_del=max_del,
        neighbor_size=neighbor_size,
        n_iters=n_iters,
        time_limit=time_limit,
        n_threads=n_threads,
    )
    tsc.io.write(df_out, f"{outfile}.gpps.CFMatrix")

    return None
