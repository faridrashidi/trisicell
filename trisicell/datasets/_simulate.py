import math
import os
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

    onconem, onconem_is_not_imported = tsc.ul.import_rpy2(
        "oncoNEM",
        "BiocManager::install('graph')\ndevtools::install_bitbucket('edith_ross/oncoNEM')\n",
    )
    if onconem_is_not_imported:
        raise RuntimeError("Unable to import a package!")

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


def splatter(noisy, ground, seed=5):
    params = tsc.ul.get_file("trisicell.datasets/data/splatter.rds")

    splatter, splatter_is_not_imported = tsc.ul.import_rpy2(
        "splatter",
        "BiocManager::install('splatter')\n",
    )
    if splatter_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    import rpy2.robjects as ro

    cmd = f"""
    suppressPackageStartupMessages({{
        library(splatter)
    }})

    seed <- strtoi({seed})
    name <- basename('{noisy}')
    groundmatrix <- read.table(file='{ground}', header=TRUE, sep='\t', row.names=1)
    counts <- colSums(groundmatrix == 1)
    # print(length(counts[counts < 2]))

    inmatrix <- read.table(file='{noisy}', header=TRUE, sep='\t', row.names=1)
    # counts <- colSums(inmatrix == 1)
    # print(counts[counts < 2])

    n <- dim(inmatrix)[1]
    m <- dim(inmatrix)[2]

    params <- readRDS(file='{params}')
    sim <- splatSimulate(params, nGenes=m, batchCells=n, seed=seed, verbose=FALSE)
    
    outmatrix <- inmatrix
    outmatrix[t(counts(sim)) < 10] <- 3

    counts <- colSums(outmatrix == 1)
    counts <- 2-counts[counts < 2]
    names <- names(counts)
    i<-0
    for (n in names) {{
        i <- i+1
        x <- which(inmatrix[,names[i]] == 1)
        if (counts[i] <= length(x)) {{
            my_sample <- sample(x, length(x), replace=FALSE)
            j <- 0
            for (s in my_sample) {{
                if (outmatrix[s,names[i]] == 3) {{
                    j <- j+1
                    outmatrix[s,names[i]] <- 1
                    if (counts[i] == j) {{
                        break
                    }}
                }}
            }}
        }}
    }}

    counts <- colSums(outmatrix == 1)
    counts <- 2-counts[counts < 2]
    names <- names(counts)
    i<-0
    for (n in names) {{
        i <- i+1
        x <- which(groundmatrix[,names[i]] == 1)
        # x2 <- which(inmatrix[,names[i]] == 0)
        # x <- intersect(x1, x2)
        if (counts[i] <= length(x)) {{
            my_sample <- sample(x, length(x), replace=FALSE)
            j <- 0
            for (s in my_sample) {{
                if (outmatrix[s,names[i]] == 0 || outmatrix[s,names[i]] == 3) {{
                    j <- j+1
                    outmatrix[s,names[i]] <- 1
                    if (counts[i] == j) {{
                        break
                    }}
                }}
            }}
        }}
    }}

    counts <- colSums(outmatrix == 1)
    # print(length(counts[counts < 2]))


    outfile <- '{noisy[:-len('.SC')]}-seed_{seed}.SC'
    write.table(x=data.frame('cellID_mutID'=rownames(outmatrix),outmatrix),
                             file=outfile, sep='\\t', quote=FALSE, row.names=FALSE)
    """
    ro.r(cmd)


def create_splatter():
    params = tsc.ul.get_file("trisicell.datasets/data/splatter.rds")

    splatter, splatter_is_not_imported = tsc.ul.import_rpy2(
        "splatter",
        "BiocManager::install('splatter')\n",
    )
    if splatter_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    import rpy2.robjects as ro

    cmd = f"""
    suppressPackageStartupMessages({{
        library(splatter)
    }})

    readcountfile <- '~/122/expr.count.tsv'
    data <- read.table(file=readcountfile, header=TRUE, sep='\t', row.names=1)
    data <- t(as.matrix(data))
    params <- splatEstimate(data)
    saveRDS(params, file='{params}')
    """
    ro.r(cmd)


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
