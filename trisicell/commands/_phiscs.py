import os

import click

import trisicell as tsc


@click.command(short_help="Run PhISCS (CSP version).")
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
def phiscsb(genotype_file, alpha, beta):
    """PhISCS-B.

    A combinatorial approach for subperfect
    tumor phylogeny reconstructionvia integrative use of
    single-cell and bulk sequencing data :cite:`PhISCS`.

    trisicell phiscsb input.SC 0.0001 0.1
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.phiscsb.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.phiscsb(df_in, alpha=alpha, beta=beta)
    tsc.io.write(df_out, f"{outfile}.phiscsb.CFMatrix")

    return None


@click.command(short_help="Run PhISCS (ILP version).")
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
def phiscsi(genotype_file, alpha, beta, time_limit, n_threads):
    """PhISCS-I.

    A combinatorial approach for subperfect
    tumor phylogeny reconstructionvia integrative use of
    single-cell and bulk sequencing data :cite:`PhISCS`.

    trisicell phiscsi input.SC 0.0001 0.1 -t 3600 -p 8
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.phiscsi.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.phiscsi(
        df_in, alpha=alpha, beta=beta, time_limit=time_limit, n_threads=n_threads
    )
    tsc.io.write(df_out, f"{outfile}.phiscsi.CFMatrix")

    return None
