import math
import random

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc


def simulate2(n_cells=10, n_muts=10, n_clones=3, alpha=0.00001, beta=0.1, missing=0):
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
    df = pd.DataFrame(dat, dtype=int)
    df.columns = [f"mut{x}" for x in df.columns]
    df.index = [f"cell{x}" for x in df.index]

    return df


def simulate(n_cells=10, n_muts=10, seed=0):
    tree = _simulate_binary_tree(n_cells, seed=seed)
    return tree


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
                    sys.exit(2)

    tsc.logg.info(f"FNs={countFN}, FPs={countFP}, NAs={countNA}")

    df_out = pd.DataFrame(data2)
    df_out.columns = df_in.columns
    df_out.index = df_in.index
    df_out.index.name = "cellIDxmutID"

    return df_out


def _helper_binary_tree(tree, root, n_leaves):
    if n_leaves == 1:
        return tree
    if n_leaves >= 2:
        tree.add_node(2 * root)
        tree.add_edge(root, 2 * root)
        tree.add_node(2 * root + 1)
        tree.add_edge(root, 2 * root + 1)
    new_n_leaves = np.random.randint(1, n_leaves)
    _helper_binary_tree(tree, 2 * root, new_n_leaves)
    _helper_binary_tree(tree, 2 * root + 1, n_leaves - new_n_leaves)


def _simulate_binary_tree(n_samples, seed=None):
    if seed:
        np.random.seed(seed)
    G = nx.DiGraph()
    _helper_binary_tree(G, 1, n_samples)
    for n in G:
        G.nodes[n]["cell"] = tsc.ul.Cell()
    return G
