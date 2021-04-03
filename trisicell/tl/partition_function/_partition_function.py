import datetime
import pickle
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import trisicell as tsc

from ._pf import count_disticnt_matrices, get_samples, get_samples_info, process_samples


def _save_samples(filename, edges_list, subtrees_list, tree_our_prob_list):
    samples_object = (edges_list, subtrees_list, tree_our_prob_list)
    with open(filename, "wb") as f:
        pickle.dump(samples_object, f)


def _load_samples(filename):
    with open(filename, "rb") as f:
        edges_list, subtrees_list, tree_our_prob_list = pickle.load(f)
    return edges_list, subtrees_list, tree_our_prob_list


def partition_function(df_input, alpha, beta, n_samples, n_batches, muts, cells):
    """Calculating the probability of a mutation seeding particular cells.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix.
    alpha : float
        False positive error rate.
    beta : float
        False negative error rate.
    n_samples : int
        Number of samples to get from the distribution (suggest: 1000)
    n_batches : int
        Number of batches to repeat the experiment (suggest: 100)
    muts : list
        The list of mutations
    cells : list
        The list of cells

    Returns
    -------
    :class:`pandas.DataFrame`
        A table of probabilities for every mutation and every batch.

    Raises
    ------
    ValueError
        Mutations must be present in the input data.
    ValueError
        Cells must be present in the input data.
    """
    df_output = pd.DataFrame(None, index=muts, columns=range(n_batches))
    s_time = time.time()
    I = df_input.values
    t1 = I * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I == 3] = 0.5

    my_muts = np.where(df_input.columns.isin(muts))[0]
    my_cells = np.where(df_input.index.isin(cells))[0]
    if len(my_muts) != len(muts):
        raise ValueError("bad muts choise!")
    if len(my_cells) != len(cells):
        raise ValueError("bad cells choise!")

    edges_list, subtrees_list, tree_our_prob_list = get_samples(P, n_samples)
    # _save_samples(
    #     './working/PF/bwes.pkl',
    #     edges_list, subtrees_list, tree_our_prob_list
    # )

    # edges_list, subtrees_list, tree_our_prob_list = _load_samples(
    #     './working/PF/bwes.pkl'
    # )

    def run(mut):
        my_mut = np.where(df_input.columns == mut)[0][0]
        pf_cond_list, tree_origin_prob_list = get_samples_info(
            P, my_cells, my_mut, n_samples, subtrees_list
        )
        estimates = process_samples(
            pf_cond_list, tree_origin_prob_list, tree_our_prob_list, n_batches
        )
        return mut, estimates

    output = Parallel(n_jobs=len(muts))(delayed(run)(mut) for mut in muts)
    for mut, estimates in output:
        df_output.loc[mut] = estimates

    # for mut in muts[:2]:
    #     mut, estimates = run(mut)
    #     df_output.loc[mut] = estimates

    e_time = time.time()
    running_time = e_time - s_time
    tsc.logg.info(f"elapsed time: {datetime.timedelta(seconds=running_time)}")

    return df_output.astype(float)
