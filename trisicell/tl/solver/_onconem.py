import os
import time
from logging import root

import numpy as np
import pandas as pd

import trisicell as tsc


def onconem(df_input, alpha, beta):
    """Running OncoNEM.

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
    from rpy2.robjects import pandas2ri

    tsc.logg.info(f"running OncoNEM with alpha={alpha}, beta={beta}")

    pandas2ri.activate()

    df_r = ro.conversion.py2rpy(df_input.replace(3, 2).T)
    ro.globalenv["df"] = df_r

    # onem = onconem.oncoNEM[new(Data=df_r, FPR=alpha, FNR=beta)
    # print(onconem.oncoNEM.rx2("new"))

    cmd = f"""
    suppressPackageStartupMessages({{
        library(oncoNEM)
        library(igraph)
    }})

    mat <- data.matrix(df)
    oNEM <- oncoNEM$new(Data=mat, FPR=as.numeric({alpha}), FNR=as.numeric({beta}))
    oNEM$search(delta=200)
    oNEM.expanded <- expandOncoNEM(oNEM, epsilon=10, delta=200,
                                   checkMax=10000, app=TRUE)
    oncoTree <-clusterOncoNEM(oNEM=oNEM.expanded, epsilon=10)
    edges <- get.edgelist(oncoTree$g)
    post <- oncoNEMposteriors(tree=oncoTree$g, clones=oncoTree$clones,
                              Data=oNEM$Data, FPR=oNEM$FPR, FNR=oNEM$FNR)
    
    post$p_mut
    """
    s_time = time.time()
    result = ro.r(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    sol_Y = np.rint(result)

    df_output = pd.DataFrame(sol_Y)
    df_output = df_output.T
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"
    df_output = df_output.astype(int)

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
