import os
import time

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc


def scite(
    df_input,
    alpha,
    beta,
    n_iters=90000,
    n_restarts=3,
    save_inter=False,
    dir_inter=".",
    base_inter=None,
    experiment=False,
):
    """Solving using SCITE.

    Tree inference for single-cell data :cite:`SCITE`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    n_iters : :obj:`int`, optional
        Number of iterations, by default 90000
    n_restarts : :obj:`int`, optional
        Number of restarts, by default 3
    save_inter : :obj:`bool`, optional
        Save input/output of the ScisTree format, by default False
    dir_inter : :obj:`str`, optional
        Directory of the output of the SCITE, by default "."
    base_inter : :obj:`str`, optional
        Basename of the output of the SCITE, by default None
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    if not experiment:
        tsc.logg.info(
            f"running SCITE with alpha={alpha}, beta={beta}, n_iters={n_iters}, "
            f"n_restarts={n_restarts}"
        )
    if not base_inter:
        tmpdir = tsc.ul.tmpdir(prefix="trisicell.", suffix=".scite", dirname=dir_inter)
    else:
        tmpdir = tsc.ul.mkdir(os.path.join(dir_inter, base_inter))

    np.savetxt(f"{tmpdir}/scite.SC.T", df_input.values.T, delimiter="\t", fmt="%1.0f")
    with open(f"{tmpdir}/scite.geneNames", "w") as fout:
        fout.write("\n".join(df_input.columns))

    scite = tsc.ul.get_file("trisicell.external/bin/scite")
    cmd = (
        f"{scite} "
        f"-i {tmpdir}/scite.SC.T "
        f"-names {tmpdir}/scite.geneNames "
        f"-n {df_input.shape[1]} "
        f"-m {df_input.shape[0]} "
        f"-ad {beta} "
        f"-fd {alpha} "
        f"-r {n_restarts} "
        "-e 0.20 "
        "-a "
        f"-l {n_iters} "
        f"-o {tmpdir}/scite > {tmpdir}/scite.log"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    with open(f"{tmpdir}/scite_ml0.gv") as fin:
        with open(f"{tmpdir}/scite_ml0_quoted.gv", "w") as fout:
            for line in fin:
                if " -> " in line:
                    line = line.strip()
                    a = line.split(" -> ")[0]
                    b = line.split(" -> ")[1].replace(";", "")
                    fout.write(f'"{a}" -> "{b}";\n')
                else:
                    fout.write(line)

    detail = {}
    with open(f"{tmpdir}/scite.log") as fin:
        for line in fin:
            line = line.strip()
            if "best value for beta:" in line:
                detail["beta"] = float(line.replace("best value for beta:", "").strip())
            if "best log score for tree:" in line:
                detail["score"] = float(
                    line.replace("best log score for tree:", "").strip()
                )

    G = nx.drawing.nx_pydot.read_dot(f"{tmpdir}/scite_ml0_quoted.gv")
    df_output = df_input.copy()
    for col in df_output.columns:
        df_output[col].values[:] = 0
    for i in range(df_output.shape[0]):
        muts = nx.shortest_path(G, source="Root", target=f"s{i}")
        muts.remove("Root")
        muts.remove(f"s{i}")
        if len(muts) > 0:
            df_output.loc[df_output.index[i], muts] = 1

    if not save_inter:
        tsc.ul.cleanup(tmpdir)

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        for k, v in detail.items():
            tsc.logg.info(f"{k}: {v}")
        return df_output
    else:
        return df_output, running_time, detail["score"], detail["beta"]


def infscite(
    df_input,
    alpha,
    beta,
    n_iters,
    n_restarts=3,
    save_inter=False,
    dir_inter=".",
    base_inter=None,
    experiment=False,
):
    tsc.logg.info(
        f"running infSCITE with alpha={alpha}, beta={beta}, n_iters={n_iters}, "
        f"n_restarts={n_restarts}"
    )
    if not base_inter:
        tmpdir = tsc.ul.tmpdir(
            prefix="trisicell_", suffix="_infscite", dirname=dir_inter
        )
    else:
        tmpdir = tsc.ul.mkdir(os.path.join(dir_inter, base_inter))

    np.savetxt(
        f"{tmpdir}/infscite.SC.T", df_input.values.T, delimiter="\t", fmt="%1.0f"
    )
    with open(f"{tmpdir}/infscite.geneNames", "w") as fout:
        fout.write("\n".join(df_input.columns))

    infscite = tsc.ul.get_file("trisicell.external/bin/infSCITE")
    cmd = (
        f"{infscite} "
        f"-i {tmpdir}/infscite.SC.T "
        f"-names {tmpdir}/infscite.geneNames "
        f"-n {df_input.shape[1]} "
        f"-m {df_input.shape[0]} "
        f"-ad {beta} "
        f"-fd {alpha} "
        f"-r {n_restarts} "
        "-z "
        "-e 0.20 "
        "-a "
        f"-l {n_iters} "
        f"-o {tmpdir}/infscite > {tmpdir}/infscite.log"
    )
    # "-rec 1 "
    # "-d "
    # "-s -p 10000 "

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    with open(f"{tmpdir}/infscite_ml0.gv") as fin:
        with open(f"{tmpdir}/infscite_ml0_quoted.gv", "w") as fout:
            for line in fin:
                if " -> " in line:
                    line = line.strip()
                    a = line.split(" -> ")[0]
                    b = line.split(" -> ")[1].replace(";", "")
                    fout.write(f'"{a}" -> "{b}";\n')
                else:
                    fout.write(line)

    detail = {}
    with open(f"{tmpdir}/infscite.log") as fin:
        for line in fin:
            line = line.strip()
            if "best value for beta:" in line:
                detail["beta"] = float(line.replace("best value for beta:", "").strip())
            if "best value for alpha:" in line:
                detail["alpha"] = float(
                    line.replace("best value for alpha:", "").strip()
                )
            if "best doublet rate:" in line:
                detail["doublet"] = float(
                    line.replace("best doublet rate:", "").strip()
                )
            if "best log score for tree:" in line:
                detail["score"] = float(
                    line.replace("best log score for tree:", "").strip()
                )

    G = nx.drawing.nx_pydot.read_dot(f"{tmpdir}/infscite_ml0_quoted.gv")
    df_output = df_input.copy()
    for col in df_output.columns:
        df_output[col].values[:] = 0
    for i in range(df_output.shape[0]):
        muts = nx.shortest_path(G, source=f"{df_input.shape[0]}", target=f"{i+1}")
        # muts.remove("Root")
        muts.remove(f"{i+1}")
        if len(muts) > 0:
            df_output.loc[df_output.index[i], muts] = 1

    if not save_inter:
        tsc.ul.cleanup(tmpdir)

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        for k, v in detail.items():
            tsc.logg.info(f"{k}: {v}")
        return df_output
    else:
        return df_output, running_time, detail["score"], detail["beta"]
