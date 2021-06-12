import os
import time

import trisicell as tsc
from trisicell.tl.solver.booster._dependencies import prepare_dependencies
from trisicell.tl.solver.booster._reconstruct_big_tree import reconstruct_big_tree
from trisicell.tl.solver.booster._subsamples import subsampling


def booster(
    df_input,
    alpha,
    beta,
    solver="SCITE",
    sample_on="muts",
    sample_size=10,
    n_samples=10,
    begin_index=0,
    n_jobs=10,
    dep_weight=50,
    time_out=120,
    n_iterations=500000,
    subsample_dir=None,
    disable_tqdm=False,
    no_subsampling=False,
    no_dependencies=False,
    no_reconstruction=False,
):
    """Trisicell-Boost solver.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        input noisy dataframe
    alpha : float
        false positive rate
    beta : float
        false negative rate
    solver : str, optional
        which tool is boosted {"SCITE", "PhISCS"}, by default "SCITE"
    sample_on : str, optional
        on which dimension is subsampled {"muts", "cells"}, by default "muts"
    sample_size : int, optional
        number of subsampled mutations or cells depends on `sample_on`, by default 10
    n_samples : int, optional
        number of samples, by default 10
    begin_index : int, optional
        start index of intermediate file names, by default 0
    n_jobs : int, optional
        number of jobs, by default 10
    dep_weight : int, optional
        weight multiplier, by default 50
    time_out : int, optional
        time out needed for PhISCS running on each instance, by default 120
    n_iterations : int, optional
        number of iterations needed for SCITE running, by default 500000
    subsample_dir : [type], optional
        for keeping the intermediate subsamples CFMatrices, by default None
    disable_tqdm : bool, optional
        disable progress bar, by default False
    no_subsampling : bool, optional
        subsampling (step 1/3) gets off, by default False
    no_dependencies : bool, optional
        dependencies calculation (step 2/3) gets off, by default False
    no_reconstruction : bool, optional
        reconstruction of big tree (step 3/3) gets off, by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        [description]

    See Also
    --------
    :func:`trisicell.tl.scite`.
    :func:`trisicell.tl.phiscsb`
    """

    if subsample_dir is not None:
        tmpdir = tsc.ul.mkdir(subsample_dir)
    else:
        # tmpdir = tsc.ul.tmpdirsys(suffix=".booster")
        # tmpdir = tmpdir.name
        tmpdir = tsc.ul.tmpdir(suffix=".booster")

    detail = {}

    s_time = time.time()
    # subsampling matrices and solving them
    if not no_subsampling:
        subsampling(
            df_input,
            alpha=alpha,
            beta=beta,
            solver=solver,
            sample_on=sample_on,
            sample_size=sample_size,
            n_samples=n_samples,
            begin_sample=begin_index,
            n_jobs=n_jobs,
            time_out=time_out,
            n_iterations=n_iterations,
            tmpdir=tmpdir,
            disable_tqdm=disable_tqdm,
        )

    # preparing dependencies file
    if not no_dependencies:
        n_muts = df_input.shape[1]
        max_num_submatrices = int(dep_weight * (n_muts ** 2) / (sample_size ** 2))
        prepare_dependencies(
            df_input.columns,
            tmpdir,
            f"{tmpdir}/_booster.dependencies",
            max_num_submatrices,
            disable_tqdm,
        )

    # building the final CFMatrix
    if not no_reconstruction:
        tsc.io.write(df_input, f"{tmpdir}/_input.SC")
        detail["TREE_SCORE"] = reconstruct_big_tree(
            f"{tmpdir}/_booster.dependencies",
            f"{tmpdir}/_input.SC",
            alpha,
            beta,
            f"{tmpdir}/_booster",
            disable_tqdm,
        )
        df_output = tsc.io.read(
            f"{tmpdir}/_booster.dnc.CFMatrix",
        )
        df_output = df_output.loc[df_input.index, df_input.columns]
    else:
        df_output = None
    e_time = time.time()
    running_time = e_time - s_time

    if subsample_dir is None:
        tsc.ul.cleanup(tmpdir)

    if df_output is not None:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        for k, v in detail.items():
            tsc.logg.info(f"{k}: {v}")

    return df_output
