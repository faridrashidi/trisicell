import os

import click
from joblib import Parallel, delayed
from tqdm import tqdm

import trisicell as tsc
from trisicell.pl._trees import _newick_info2_mutation_list


def run_scistree(df_in, alpha, beta, outfile):
    tsc.settings.logfile = f"{outfile}/fn_{beta}-fp_{alpha}.log"
    df_out = tsc.tl.scistree(df_in, alpha, beta)
    tsc.io.write(df_out, f"{outfile}/fn_{beta}-fp_{alpha}.CFMatrix")

    tree = tsc.ul.to_tree(df_out)
    newick, info2, mutations = _newick_info2_mutation_list(tree)
    with open(f"{outfile}/fn_{beta}-fp_{alpha}.newick", "w") as fout:
        fout.write(newick + "\n")
    info2.to_csv(f"{outfile}/fn_{beta}-fp_{alpha}.info2", index=None)


@click.command(short_help="Grid search for all parameters.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def search(genotype_file):
    """Grid search for all parameters of alpha and beta.

    trisicell search input.SC
    """

    tsc.settings.verbosity = "info"

    outfile = os.path.splitext(genotype_file)[0]
    tsc.ul.mkdir(outfile)

    df_in = tsc.io.read(genotype_file)

    betas = [0.1, 0.2, 0.3, 0.4]
    alphas = [0.1, 0.01, 0.001, 0.0001, 0.00001]
    n_samples = len(betas) * len(alphas)

    with tsc.ul.tqdm_joblib(
        tqdm(
            ascii=True,
            ncols=100,
            desc="GRID SEARCHING:",
            total=n_samples,
            position=0,
        )
    ) as progress_bar:
        output = Parallel(n_jobs=n_samples)(
            delayed(run_scistree)(df_in, alpha, beta, outfile)
            for alpha in alphas
            for beta in betas
        )

    return None
