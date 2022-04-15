import os

import click

import trisicell as tsc


@click.command(short_help="Run PhISCS-BnB.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.option(
    "--bounding",
    "-b",
    default="real",
    type=click.Choice(["real", "simulated"]),
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--time_limit",
    "-t",
    default=None,
    type=int,
    show_default=True,
    help="Time limit.",
)
def bnb(genotype_file, bounding, time_limit):
    """PhISCS-BnB.

    A fast branch and bound algorithm for
    the perfect tumor phylogeny reconstruction problem :cite:`PhISCS-BnB`.

    trisicell bnb input.SC -b real
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"
    tsc.settings.logfile = f"{outfile}.bnb.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.bnb(df_in, bounding=bounding, time_limit=time_limit)
    tsc.io.write(df_out, f"{outfile}.bnb.CFMatrix")

    return None
