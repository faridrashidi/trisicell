import datetime
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import trisicell as tsc
from trisicell.tl.partition_function._pf import (
    get_samples,
    get_samples_info,
    process_samples,
)


def partition_function(df_input, alpha, beta, n_samples, n_batches, muts, cells):
    """Calculate the probability of a mutation seeding particular cells.

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
    """

    df_output = pd.DataFrame(None, index=muts, columns=range(n_batches))
    s_time = time.time()
    I_mtr = df_input.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    my_muts = np.where(df_input.columns.isin(muts))[0]
    my_cells = np.where(df_input.index.isin(cells))[0]
    if len(my_muts) != len(muts):
        tsc.logg.error("bad muts choise!")
    if len(my_cells) != len(cells):
        tsc.logg.error("bad cells choise!")

    _, subtrees_list, tree_our_prob_list = get_samples(P, n_samples)

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

    e_time = time.time()
    running_time = e_time - s_time
    tsc.logg.info(f"elapsed time: {datetime.timedelta(seconds=running_time)}")

    return df_output.astype(float)
