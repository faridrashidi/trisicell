import math
import random

import numpy as np
import pandas as pd

import trisicell as tsc


def simulate(n_cells=10, n_muts=10, n_clones=3, alpha=0.00001, beta=0.1, missing=0):
    """Simulate single-cell noisy genotype matrix.

    This function is using :cite:`OncoNEM`.

    Parameters
    ----------
    n_cells : :obj:`int`, optional
        Number of cells, by default 10
    n_muts : :obj:`int`, optional
        Number of mutations, by default 10
    n_clones : :obj:`int`, optional
        Number of clones, by default 3
    alpha : :obj:`float`, optional
        False positive rate, by default 0.00001
    beta : :obj:`float`, optional
        False negative rate, by default 0.1
    missing : :obj:`int`, optional
        Missing entry rate, by default 0

    Returns
    -------
    :class:`pandas.DataFrame`
        A genotype matrix where 0 is absent, 1 is present and 3 is missing.
    """

    # TODO: replace

    onconem, onconem_is_not_imported = tsc.ul.import_rpy2(
        "oncoNEM",
        "BiocManager::install('graph')\ndevtools::install_bitbucket('edith_ross/oncoNEM')\n",
    )
    if onconem_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    dat = onconem.simulateData(
        N_cells=n_cells,
        N_normalContam=0,
        N_clones=n_clones,
        N_unobs=0,
        N_sites=n_muts,
        FPR=alpha,
        FNR=beta,
        p_missing=missing,
        randomizeOrder=False,
    )

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        dat = ro.conversion.rpy2py(dat.rx2("D"))
    dat[dat == 2] = 3
    df = pd.DataFrame(dat.T, dtype=int)
    df.columns = [f"mut{x}" for x in df.columns]
    df.index = [f"cell{x}" for x in df.index]

    return df


def add_noise(df_in, alpha, beta, missing):
    """Add noise to the input genotype matrix.

    These noise includes:
    1) False positive errors (alpha)
    2) False negative errors (beta)
    3) Missing entry errors (missing)

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    missing : :obj:`float`
        Missing entry error rate.

    Returns
    -------
    :class:`pandas.DataFrame`
        A noisy genotype matrix where 0 is absent, 1 is present and 3 is missing.
    """

    def toss(p):
        return True if random.random() < p else False

    data = df_in.values
    n, m = df_in.shape
    data2 = -1 * np.ones(shape=(n, m)).astype(int)
    countFP = 0
    countFN = 0
    countNA = 0
    countOneZero = 0
    indexNA = []
    changedBefore = []
    for i in range(n):
        for j in range(m):
            indexNA.append([i, j])
            countOneZero = countOneZero + 1
    random.shuffle(indexNA)
    nas = math.ceil(countOneZero * missing)
    for i in range(int(nas)):
        [a, b] = indexNA[i]
        changedBefore.append([a, b])
        data2[a][b] = 3
        countNA = countNA + 1
    for i in range(n):
        for j in range(m):
            if data2[i][j] != 3:
                if data[i][j] == 1:
                    if toss(beta):
                        data2[i][j] = 0
                        countFN = countFN + 1
                    else:
                        data2[i][j] = data[i][j]
                elif data[i][j] == 0:
                    if toss(alpha):
                        data2[i][j] = 1
                        countFP = countFP + 1
                    else:
                        data2[i][j] = data[i][j]
                else:
                    tsc.logg.error("Wrong Input")

    # tsc.logg.info(f"FNs={countFN}, FPs={countFP}, NAs={countNA}")

    df_out = pd.DataFrame(data2)
    df_out.columns = df_in.columns
    df_out.index = df_in.index
    df_out.index.name = "cellIDxmutID"

    return df_out


def add_doublets(df_ground, df_noisy, alpha, beta, missing, doublet):
    df_doublet = df_noisy.copy()
    doublet_cells = []
    for _ in range(int(doublet * df_ground.shape[0])):
        r1 = np.random.choice(df_ground.index, replace=False, size=1)
        while r1 in doublet_cells:
            r1 = np.random.choice(df_ground.index, replace=False, size=1)
        doublet_cells.append(r1)
        r2 = np.random.choice(df_ground.index, replace=False, size=1)
        df_doublet.loc[r1] = 1 * np.logical_or(df_ground.loc[r1], df_ground.loc[r2])
        df_doublet.loc[r1] = tsc.datasets.add_noise(
            df_doublet.loc[r1], alpha=alpha, beta=beta, missing=missing
        )
    return df_doublet
