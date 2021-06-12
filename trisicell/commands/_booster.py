import os

import click

import trisicell as tsc


@click.command(short_help="Run Booster.")
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
    "--solver",
    default="scite",
    type=click.Choice(["scite", "phiscs", "scistree"]),
    show_default=True,
    help="Solver of the booster.",
)
@click.option(
    "--sample_on",
    default="muts",
    type=click.Choice(["muts", "cells"]),
    show_default=True,
    help="Sampling on `muts` or `cells`.",
)
@click.option(
    "--sample_size",
    default=10,
    type=int,
    show_default=True,
    help="The size of samples i.e. the number of muts or cells.",
)
@click.option(
    "--n_samples",
    default=10,
    type=int,
    show_default=True,
    help="The number of subsamples.",
)
@click.option(
    "--begin_index",
    default=0,
    type=int,
    show_default=True,
    help="ID of the start subsample name.",
)
@click.option(
    "--n_jobs",
    default=0,
    type=int,
    show_default=True,
    help="Number of parallel jobs to do subsampleing.",
)
@click.option(
    "--dep_weight",
    default=50,
    type=int,
    show_default=True,
    help="""Weight for how many subsamples to be used in dependencies calculation.""",
)
@click.option(
    "--time_out",
    default=120,
    type=int,
    show_default=True,
    help="Timeout of solving allowance.",
)
@click.option(
    "--n_iterations",
    default=500000,
    type=int,
    show_default=True,
    help="SCITE number of iterations.",
)
@click.option(
    "--subsample_dir",
    default=None,
    type=str,
    show_default=True,
    help="A path to the subsamples directory.",
)
@click.option(
    "--disable_tqdm",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Disable showing the tqdm progress.",
)
@click.option(
    "--no_subsampling",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="No running of subsampling (step 1/3).",
)
@click.option(
    "--no_dependencies",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="No running of subsampling (step 2/3).",
)
@click.option(
    "--no_reconstruction",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="No running of reconstruction (step 3/3).",
)
def booster(
    genotype_file,
    alpha,
    beta,
    solver,
    sample_on,
    sample_size,
    n_samples,
    begin_index,
    n_jobs,
    dep_weight,
    time_out,
    n_iterations,
    subsample_dir,
    disable_tqdm,
    no_subsampling,
    no_dependencies,
    no_reconstruction,
):
    """Divide and Conquer

    For doing all 3 steps:

    trisicell booster input.SC 0.001 0.1
    --solver scite --n_samples 200 --sample_size 15
    --n_jobs 4 --n_iterations 10000 --dep_weight 50
    --subsample_dir . --begin_index 0

    For doing only the last step:

    trisicell booster input.SC 0.001 0.1
    --dep_weight 50 --subsample_dir PATH_TO_SUBSAMPLES_FOLDER
    --no_subsampling --no_dependencies
    """

    dirbase = tsc.ul.dirbase(genotype_file)
    # dirname, basename = tsc.ul.dir_base(genotype_file)

    tsc.settings.verbosity = "info"
    if not no_reconstruction:
        tsc.settings.logfile = f"{dirbase}.booster.log"

    df_in = tsc.io.read(genotype_file)
    df_out = tsc.tl.booster(
        df_in,
        alpha=alpha,
        beta=beta,
        solver=solver,
        sample_on=sample_on,
        sample_size=sample_size,
        n_samples=n_samples,
        begin_index=begin_index,
        n_jobs=n_jobs,
        dep_weight=dep_weight,
        time_out=time_out,
        n_iterations=n_iterations,
        subsample_dir=subsample_dir,
        disable_tqdm=disable_tqdm,
        no_subsampling=no_subsampling,
        no_dependencies=no_dependencies,
        no_reconstruction=no_reconstruction,
    )
    if not no_reconstruction:
        tsc.io.write(df_out, f"{dirbase}.booster.CFMatrix")

    return None
