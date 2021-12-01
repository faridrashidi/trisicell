import time

import pandas as pd

import trisicell as tsc
from trisicell.external.gpps import gpps_hc, gpps_ilp


def gpps(
    df_input,
    alpha,
    beta,
    k_dollo=0,
    max_del=-1,
    neighbor_size=30,
    n_iters=100,
    time_limit=86400,
    n_threads=1,
):
    """Solving using gpps.

    an ILP-based approach for inferring cancer progression with mutation losses from
    single cell data :cite:`gpps`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    k_dollo : :obj:`int`, optional
        k for Dollo model, by default 0
    max_del : :obj:`int`, optional
        Maximum number of deletion allowed, by default -1
    neighbor_size : :obj:`int`, optional
        Hill climbing neighborhood size, by default 30
    n_iters : :obj:`int`, optional
        Hill climbing maximum iterations, by default 100
    time_limit : :obj:`int`, optional
        Time limit (in seconds), by default 86400
    n_threads : :obj:`int`, optional
        Number of threads, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    tsc.logg.info(
        f"running gpps with alpha={alpha}, beta={beta}, k_dollo={k_dollo}, "
        f"max_del={max_del}, neighbor_size={neighbor_size}, n_iters={n_iters}, "
        f"time_limit={time_limit}, n_threads={n_threads}"
    )

    cells = list(df_input.index)
    snvs = list(df_input.columns)

    s_time = time.time()
    ilp_matrix = gpps_ilp(
        df_input.values,
        alpha=beta,  # gpps takes a as false-negative and b as false-positive
        beta=alpha,
        k_dollo=k_dollo,
        max_del=max_del,
        time_limit=time_limit,
        n_threads=n_threads,
    )
    ilp_matrix = pd.DataFrame(ilp_matrix)

    output_matrix = gpps_hc(
        df_input.values,
        ilp_matrix,
        alpha=beta,
        beta=alpha,
        k_dollo=k_dollo,
        mut_names=snvs,
        ns=neighbor_size,
        mi=n_iters,
    )
    e_time = time.time()
    running_time = e_time - s_time

    df_output = pd.DataFrame(output_matrix)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
