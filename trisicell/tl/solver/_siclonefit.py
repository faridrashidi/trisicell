import glob
import os
import time

import numpy as np
import pandas as pd

import trisicell as tsc


def siclonefit(
    df_input,
    alpha,
    beta,
    n_restarts=3,
    n_iters=500,
    burnin=100,
    return_tree=False,
    experiment=False,
):
    """Solving using SiCloneFit.

    Bayesian inference of population structure, genotype, and phylogeny of tumor clones
    from single-cell genome sequencing data :cite:`SiCloneFit`.

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
    n_restarts : :obj:`int`, optional
        Number of restarts, by default 3
    n_iters : :obj:`int`, optional
        Number of iterations, by default 90000
    return_tree : :obj:`bool`, optional
        Return the inferred cell-lineage tree, by default False
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    executable = tsc.ul.executable("SiCloneFiTComplete.jar", "SiCloneFit")

    if not experiment:
        tsc.logg.info(
            f"running SiCloneFit with alpha={alpha}, beta={beta}, n_iters={n_iters}"
        )

    tmpdir = tsc.ul.tmpdirsys(suffix=".siclonefit")

    df_input.T.reset_index(drop=True).to_csv(
        f"{tmpdir.name}/siclonefit.input", sep=" ", header=None
    )
    with open(f"{tmpdir.name}/siclonefit.cellnames", "w") as fout:
        fout.write(" ".join(df_input.index))
    with open(f"{tmpdir.name}/siclonefit.genenames", "w") as fout:
        fout.write(" ".join(df_input.columns))
    I_mtr = df_input.values

    cmd = (
        f"java -jar {executable} "
        f"-m {df_input.shape[0]} "
        f"-n {df_input.shape[1]} "
        f"-ipMat {tmpdir.name}/siclonefit.input "
        f"-fp {alpha} "
        f"-fn {beta} "
        "-df 0 "
        f"-missing {np.sum(I_mtr == 3)/(I_mtr.size)} "
        f"-iter {n_iters} "
        f"-cellNames {tmpdir.name}/siclonefit.cellnames "
        f"-geneNames {tmpdir.name}/siclonefit.genenames "
        f"-r {n_restarts} "
        f"-burnin {burnin} "
        # "-recurProb 0 "
        # "-delProb 0 "
        # "-LOHProb 0 "
        # "-doublet "
        # "-printIter "
        # "-treeIter "
        f"-outDir {tmpdir.name} > {tmpdir.name}/siclonefit.log"
    )
    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    out_dir = glob.glob(f"{tmpdir.name}/*samples/best")[0]

    df_output = pd.read_csv(
        f"{out_dir}/best_MAP_predicted_genotype.txt",
        sep=" ",
        header=None,
        index_col=0,
    ).T
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    with open(f"{out_dir}/best_MAP_tree.txt") as fin:
        tree = fin.readline().strip()

    tmpdir.cleanup()

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        if return_tree:
            return df_output, tree
        else:
            return df_output
    else:
        is_cf = tsc.ul.is_conflict_free_gusfield(df_output)
        nll = tsc.ul.calc_nll_matrix(df_input, df_output, alpha, beta)
        return df_output, running_time, is_cf, nll
