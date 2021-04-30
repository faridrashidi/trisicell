import numba
import numpy as np
import pandas as pd
import scipy as sp
from scipy.cluster.hierarchy import cut_tree, dendrogram, fcluster, linkage
from sklearn.metrics import pairwise_distances

import trisicell as tsc
from trisicell.external._betabinom import pmf_BetaBinomial


@numba.jit(nopython=True)
def _l1_ignore_na(a, b):
    a[a == 3] = np.nan
    b[b == 3] = np.nan
    return np.nanmean(np.abs(a - b))


def dist_l1_ignore_na(I, n_jobs=1):
    dist = pairwise_distances(
        I, metric=_l1_ignore_na, force_all_finite="allow-nan", n_jobs=n_jobs
    )
    np.fill_diagonal(dist, 0)
    return dist


# https://gist.github.com/FedericoV/0e7d6d8c8794a99a7a42
@numba.jit(nopython=True)
def _cosine_ignore_na(u, v):
    m = u.shape[0]
    udotv = 0
    u_norm = 0
    v_norm = 0
    for i in range(m):
        if (np.isnan(u[i])) or (np.isnan(v[i])):
            continue
        udotv += u[i] * v[i]
        u_norm += u[i] * u[i]
        v_norm += v[i] * v[i]
    u_norm = np.sqrt(u_norm)
    v_norm = np.sqrt(v_norm)
    if (u_norm == 0) or (v_norm == 0):
        ratio = 1.0
    else:
        ratio = 1 - udotv / (u_norm * v_norm)
    if ratio < 0:
        return 0
    return ratio


def dist_cosine_ignore_na(I, n_jobs=1):
    dist = pairwise_distances(
        I, metric=_cosine_ignore_na, force_all_finite="allow-nan", n_jobs=n_jobs
    )
    np.fill_diagonal(dist, 0)
    return dist


def _dist_dendro(T, V, I):
    PROB_SEQ_ERROR = 0.001

    def logSum_1(x, y):
        big = np.copy(x)
        big[x < y] = y[x < y]
        small = np.copy(x)
        small[x >= y] = y[x >= y]
        tmp = big + np.log(1 + np.exp(small - big))
        # tmp[np.bitwise_and(x==-np.inf, y==-np.inf)] = -np.inf
        # tmp[np.bitwise_and(x==np.inf, y==np.inf)] = np.inf
        return tmp

    D = np.divide(V, T)

    Mu = np.nanmean(D, axis=0)
    Var = np.nanvar(D, axis=0, ddof=1)
    a = ((1 - Mu) * Mu / Var - 1) * Mu
    b = ((1 - Mu) * Mu / Var - 1) * (1 - Mu)
    bad_muts = (
        (a <= 0) | (b <= 0) | np.isnan(a) | np.isnan(b) | np.isinf(a) | np.isinf(b)
    )
    V = V[:, ~bad_muts]
    T = T[:, ~bad_muts]
    D = D[:, ~bad_muts]
    I = I[:, ~bad_muts]
    a = a[~bad_muts]
    b = b[~bad_muts]

    lPz0 = np.zeros(T.shape, dtype=np.float64)
    lPz1 = np.zeros(T.shape, dtype=np.float64)
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            if T[i, j] != 0:
                lPz0[i, j] = np.log(sp.stat.binom.pmf(V[i, j], T[i, j], PROB_SEQ_ERROR))
                lPz1[i, j] = np.log(pmf_BetaBinomial(V[i, j], T[i, j], a[j], b[j]))

    Pg = np.sum(I == 1, axis=0) / I.shape[0]
    lPg = np.log(Pg)
    l1Pg = np.log(1 - Pg)
    lupiall = logSum_1(lPz0 + l1Pg, lPz1 + lPg)

    dist = np.zeros((T.shape[0], T.shape[0]), dtype=np.float64)
    for i in range(T.shape[0]):
        ldowni = logSum_1(lPz0[i, :] + lPz0 + l1Pg, lPz1[i, :] + lPz1 + lPg)
        lupi = logSum_1(lupiall[i, :] + lupiall, ldowni)
        dist[i, :] = np.sum(lupi - ldowni, axis=1)

    dist = dist - np.min(dist) + 1
    return dist, bad_muts


def dist_dendro(adata):
    T = adata.layers["total"]
    V = adata.layers["mutant"]
    G = adata.layers["genotype"]
    G[(G == 1) | (G == 3)] = 1
    G[G == 2] = 0
    dist, bad_muts = _dist_dendro(T, V, G)
    tsc.pp.remove_mut_by_list(adata, bad_muts)
    tsc.logg.info(f"{sum(bad_muts)} mutations filtered")
    return dist


def hclustering(df, metric="l1", method="ward"):
    """Hierarchical clustering.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        The genotype matrix.
    metric: :obj:`str`, optional
        The metric option. Can be:

            - `l1`
            - `cosine`
    method : :obj:`str`, optional
        The method for the hierarchical clustering, by default "ward"

    Returns
    -------
    :obj:`dict`
        A dictionary in which keys are the number of clusters and
        values are the cluster labels for each item.
    """

    if metric == "l1":
        dist = dist_l1_ignore_na(df.values)
    elif metric == "cosine":
        dist = dist_cosine_ignore_na(df.values)
    else:
        raise ValueError("Wroing `metric` choice!")
    clust = linkage(dist[np.triu_indices(dist.shape[0], 1)], method=method)
    clusters = {}

    for i in range(2, dist.shape[0]):
        fc = fcluster(clust, i, criterion="maxclust")
        clusters[i] = pd.Series(fc, index=df.index)

    return clusters
