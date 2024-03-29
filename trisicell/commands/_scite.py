import os

import click
import numpy as np
from joblib import Parallel, delayed

import trisicell as tsc


@click.command(short_help="Run SCITE.")
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
    default=1000000,
    type=int,
    show_default=True,
    help="Number of iterations.",
)
@click.option(
    "--n_restarts",
    "-r",
    default=3,
    type=int,
    show_default=True,
    help="Number of restarts.",
)
@click.option(
    "--experiment",
    "-e",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Is in experiment mode.",
)
@click.option(
    "--time_limit",
    "-t",
    default=86400,
    type=float,
    show_default=True,
    help="Time limit for the experiment part (in seconds).",
)
@click.option(
    "--smooth_rate",
    "-s",
    default=2,
    type=float,
    show_default=True,
    help="Smooth rate for the experiment part.",
)
def scite(
    genotype_file, alpha, beta, n_iters, n_restarts, experiment, time_limit, smooth_rate
):
    """SCITE.

    Tree inference for single-cell data :cite:`SCITE`.

    trisicell scite input.SC 0.0001 0.1 -l 1000000 -r 3 -e -t 86400 -s 2
    """

    outfile = os.path.splitext(genotype_file)[0]

    tsc.settings.verbosity = "info"

    df_in = tsc.io.read(genotype_file)
    if not experiment:
        tsc.settings.logfile = f"{outfile}.scite.log"
        df_out = tsc.tl.scite(
            df_in,
            alpha=alpha,
            beta=beta,
            n_iters=n_iters,
            n_restarts=n_restarts,
        )
        tsc.io.write(df_out, f"{outfile}.scite.CFMatrix")
    else:
        tsc.settings.logfile = f"{outfile}.scite.log"
        df_out, running_time, _, _ = tsc.tl.scite(
            df_in,
            alpha=alpha,
            beta=beta,
            n_iters=30000,
            n_restarts=1,
            experiment=True,
        )
        n_iters = int(smooth_rate * 30000 * time_limit / running_time)

        def run(i):
            do, r, s, b = tsc.tl.scite(
                df_in,
                alpha=alpha,
                beta=beta,
                n_iters=n_iters,
                n_restarts=1,
                experiment=True,
            )
            return do, r, s, b

        output = Parallel(n_jobs=n_restarts)(delayed(run)(i) for i in range(n_restarts))

        scores = [x[2] for x in output]
        betas = [x[3] for x in output]
        best_i = np.argmax(scores)
        df_out = output[best_i][0]

        tsc.ul.stat(df_in, df_out, alpha, beta, output[best_i][1])
        tsc.logg.info(f"score: {output[best_i][2]}")
        tsc.logg.info(f"beta: {output[best_i][3]}")
        tsc.logg.info(f"n_iters: {n_iters}")
        tsc.logg.info(f"scores: {','.join(list(map(str, scores)))}")
        tsc.logg.info(f"betas: {','.join(list(map(str, betas)))}")
        tsc.logg.info(f"picked: {best_i}")

        tsc.io.write(df_out, f"{outfile}.scite.CFMatrix")

    return None
