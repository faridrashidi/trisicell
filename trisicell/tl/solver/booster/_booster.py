import os

import pandas as pd

import trisicell as tsc
from trisicell.tl.solver.booster._dependencies import prepare_dependencies
from trisicell.tl.solver.booster._subsamples import subsampling


def booster(
    df_input,
    alpha,
    beta,
    solver="PhISCS",
    sample_on="muts",
    sample_size=10,
    n_samples=10,
    begin_sample=0,
    n_jobs=10,
    time_out=120,
    save_inter=True,
    dir_inter=".",
    base_inter=None,
    disable_tqdm=False,
    weight=50,
    no_subsampling=False,
    no_dependencies=False,
):
    """Divide and Conquer Booster solver.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        input noisy dataframe
    alpha : float
        false positive rate
    beta : float
        false negative rate
    solver : :obj:`str`, optional
        [description], by default "PhISCS"
    sample_size : :obj:`int`, optional
        [description], by default 10
    n_samples : :obj:`int`, optional
        [description], by default 10
    begin_sample : :obj:`int`, optional
        [description], by default 0
    n_jobs : :obj:`int`, optional
        [description], by default 10
    time_out : :obj:`int`, optional
        [description], by default 120
    save_inter : :obj:`bool`, optional
        [description], by default True
    dir_inter : :obj:`str`, optional
        [description], by default "."
    base_inter : :obj:`str`, optional
        [description], by default None
    disable_tqdm : :obj:`bool`, optional
        [description], by default False
    weight : :obj:`int`, optional
        [description], by default 50

    Returns
    -------
    :class:`pandas.DataFrame`
        [description]


    See Also
    --------
    :func:`trisicell.tl.scite`.
    """

    if not base_inter:
        tmpdir = tsc.ul.tmpdir(suffix=".booster", dirname=dir_inter)
    else:
        tmpdir = tsc.ul.mkdir(os.path.join(dir_inter, base_inter))

    #### subsampling matrices and solving them
    if not no_subsampling:
        subsampling(
            df_input,
            alpha=alpha,
            beta=beta,
            solver=solver,
            sample_on=sample_on,
            sample_size=sample_size,
            n_samples=n_samples,
            begin_sample=begin_sample,
            n_jobs=n_jobs,
            time_out=time_out,
            tmpdir=tmpdir,
            disable_tqdm=disable_tqdm,
        )

    #### preparing dependencies file
    if not no_dependencies:
        n_muts = df_input.shape[1]
        max_num_submatrices = int(weight * (n_muts ** 2) / (sample_size ** 2))
        prepare_dependencies(
            df_input.columns,
            tmpdir,
            f"{tmpdir}/_booster.dependencies",
            max_num_submatrices,
            disable_tqdm,
        )

    #### building the final cfmatrix

    df_output = pd.DataFrame(df_input.values)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    if not save_inter:
        tsc.ul.cleanup(tmpdir)

    # tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
