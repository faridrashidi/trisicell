import pandas as pd

import trisicell as tsc


def huntress(df_input_filepath, alpha, beta, kind, n_threads=1, save_inter=False):
    tsc.logg.info(
        f"running HUNTRESS with alpha={alpha}, beta={beta}, kind={kind}, n_threads={n_threads}"
    )
    tmpdir = tsc.ul.tmpdir(prefix="trisicell.", suffix=".huntress", dirname=".")

    # TODO remove tmpdir and directly work with python api.
    from ._huntress_21_03_19 import Reconstruct

    running_time = 0
    if kind == "dna":
        fn_conorm = 0.1
        fp_conorm = fn_conorm * alpha / beta
        fnfp_conorm = fn_conorm / fp_conorm

        running_time = Reconstruct(
            df_input_filepath,
            f"{tmpdir}/huntress",
            n_proc=n_threads,
            fnfp=51,
            post_fn=fn_conorm,
            post_fp=fp_conorm,
        )
        df_output = pd.read_table(f"{tmpdir}/huntress_optH_.CFMatrix", index_col=0)
    elif kind == "fn":
        running_time = Reconstruct(
            df_input_filepath,
            f"{tmpdir}/huntress.CFMatrix",
            Algchoice="FN",
            n_proc=n_threads,
        )
        df_output = pd.read_table(f"{tmpdir}/huntress.CFMatrix", index_col=0)
    elif kind == "rna":
        running_time = Reconstruct(
            df_input_filepath, f"{tmpdir}/huntress", postprocessing=1, n_proc=n_threads
        )
        df_output = pd.read_table(
            f"{tmpdir}/huntress_optH_.processed.CFMatrix", index_col=0
        )

    df_input = pd.read_table(df_input_filepath, index_col=0)
    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    if not save_inter:
        tsc.ul.cleanup(tmpdir)
    return df_output
