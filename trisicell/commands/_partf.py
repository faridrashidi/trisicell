import os
import pickle
import random
import string

import click
import numpy as np
from joblib import Parallel, delayed

import trisicell as tsc
from trisicell.tl.partition_function._pf import draw_sample_clt


@click.command(short_help="Get samples or calculate for PartF.")
# @click.argument(
#     "kind",
#     required=True,
#     type=click.Choice(["samples", "calculate"]),
# )
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
    "--n_samples",
    default=1000,
    type=int,
    show_default=True,
    help="Number of samples.",
)
@click.option(
    "--n_threads",
    default=8,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def partf(genotype_file, alpha, beta, n_samples, n_threads):
    """Get samples for calculating the partition function.

    trisicell partf input.SC 0.0001 0.1
    """

    tsc.settings.verbosity = "info"

    df_input = tsc.io.read(genotype_file)
    outfile = os.path.splitext(genotype_file)[0]
    tsc.ul.mkdir(f"{outfile}.partf.samples")

    I_mtr = df_input.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5
    P = P.astype(np.float128)

    edges_list = []
    subtrees_list = []
    tree_our_prob_list = []

    def run():
        return draw_sample_clt(P, False, c=1, coef=10)

    output = Parallel(n_jobs=n_threads)(delayed(run)() for i in range(0, n_samples))

    for edges, subtrees, prior_prob in output:
        edges_list.append(edges)
        subtrees_list.append(subtrees)
        tree_our_prob_list.append(prior_prob)
    samples_object = (edges_list, subtrees_list, tree_our_prob_list)
    filename = "".join(random.choice(string.ascii_letters) for _ in range(16))
    filename = f"{outfile}.partf.samples/{filename}.pkl"
    with open(filename, "wb") as f:
        pickle.dump(samples_object, f)

    # for filename in glob.glob(f"{outfile}.partf.samples/*.pkl"):
    #     with open(filename, "rb") as f:
    #         edges_l, subtrees_l, tree_our_prob_l = pickle.load(f)
    #         edges_list += edges_l
    #         subtrees_list += subtrees_l
    #         tree_our_prob_list += tree_our_prob_l

    return None
