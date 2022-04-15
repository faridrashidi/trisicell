import pandas as pd

import trisicell as tsc
from trisicell.tl.solver.huntress._huntress import Reconstruct


def huntress(df_input_filepath, alpha, beta, n_threads=1):
    """Solving using HUNTRESS.

    HUNTRESS: Provably fast intratumor heterogeneity inference from single-cell
    sequencing data

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
    n_threads : :obj:`int`, optional
        Number of threads, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    tsc.logg.info(
        f"running HUNTRESS with alpha={alpha}, beta={beta}, n_threads={n_threads}"
    )

    tmpdir = tsc.ul.tmpdirsys(suffix=".huntress")

    running_time = 0
    if alpha == 0:
        running_time = Reconstruct(
            df_input_filepath,
            f"{tmpdir.name}/huntress.CFMatrix",
            Algchoice="FN",
            n_proc=n_threads,
        )
        df_output = pd.read_table(f"{tmpdir.name}/huntress.CFMatrix", index_col=0)
    else:
        fn_conorm = 0.1
        fp_conorm = fn_conorm * alpha / beta
        running_time = Reconstruct(
            df_input_filepath,
            f"{tmpdir.name}/huntress",
            n_proc=n_threads,
            fnfp=51,
            post_fn=fn_conorm,
            post_fp=fp_conorm,
        )
        df_output = pd.read_table(f"{tmpdir.name}/huntress.CFMatrix", index_col=0)

    df_input = pd.read_table(df_input_filepath, index_col=0)
    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    tmpdir.cleanup()

    return df_output
