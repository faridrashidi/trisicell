import os
import time
from logging import root

import numpy as np
import pandas as pd

import trisicell as tsc


def onconem(df_input, alpha, beta):
    """Solving using OncoNEM.

    Inferring tumor evolution from single-cell sequencing data :cite:`OncoNEM`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix where 0 is absent, 1 is present and 3 is missing.
    alpha : :obj:`float`
        False positive rate.
    beta : :obj:`float`
        False negative rate.

    Returns
    -------
    :class:`pandas.DataFrame`
        Conflict-free genotype matrix where 0 is absent, 1 is present.
    """

    onconem, onconem_is_not_imported = tsc.ul.import_rpy2(
        "oncoNEM",
        "BiocManager::install('graph')\ndevtools::install_bitbucket('edith_ross/oncoNEM')\n",
    )
    if onconem_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr

    igraph = importr("igraph")

    tsc.logg.info(f"running OncoNEM with alpha={alpha}, beta={beta}")

    with ro.conversion.localconverter(ro.default_converter + numpy2ri.converter):
        data = ro.conversion.py2rpy(df_input.replace(3, 2).T.values)

    s_time = time.time()
    onem = onconem.oncoNEM(Data=data, FPR=alpha, FNR=beta)
    ro.globalenv["onem"] = onem
    ro.r("onem$search(delta=200)")
    onem_expanded = onconem.expandOncoNEM(
        onem, epsilon=10, delta=200, checkMax=10000, app=True
    )
    onco_tree = onconem.clusterOncoNEM(oNEM=onem_expanded, epsilon=10)
    edges = igraph.get_edgelist(onco_tree.rx2["g"])
    post = onconem.oncoNEMposteriors(
        tree=onco_tree.rx2["g"],
        clones=onco_tree.rx2["clones"],
        Data=data,
        FPR=alpha,
        FNR=beta,
    )
    sol_Y = np.rint(post.rx2["p_mut"])
    e_time = time.time()
    running_time = e_time - s_time

    df_output = pd.DataFrame(sol_Y)
    df_output = df_output.T
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"
    df_output = df_output.astype(int)

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
