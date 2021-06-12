import os
import time

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc
from trisicell.external._scite import run_scite


def scite(
    df_input,
    alpha,
    beta,
    n_iters=90000,
    n_restarts=3,
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

    tmpdir = tsc.ul.tmpdirsys(suffix=".scite")

    np.savetxt(
        f"{tmpdir.name}/scite.SC.T", df_input.values.T, delimiter="\t", fmt="%1.0f"
    )
    with open(f"{tmpdir.name}/scite.geneNames", "w") as fout:
        fout.write("\n".join(df_input.columns))

    """
    tmpdir = "/data/frashidi/test"
    cmd = [
        "scite",
        "-i",
        f"{tmpdir.name}/scite.SC.T",
        "-names",
        f"{tmpdir.name}/scite.geneNames",
        "-n",
        f"{df_input.shape[1]}",
        "-m",
        f"{df_input.shape[0]}",
        "-ad",
        f"{beta}",
        "-fd",
        f"{alpha}",
        "-r",
        f"{n_restarts}",
        "-e",
        "0.20",
        "-a",
        "-l",
        f"{n_iters}",
        "-o",
        f"{tmpdir.name}/scite.output",
    ]
    print(" ".join(cmd))
    run_scite(cmd)
    """

    scite = tsc.ul.get_file("trisicell.external/bin/scite")
    cmd = (
        f"{scite} "
        f"-i {tmpdir.name}/scite.SC.T "
        f"-names {tmpdir.name}/scite.geneNames "
        f"-n {df_input.shape[1]} "
        f"-m {df_input.shape[0]} "
        f"-ad {beta} "
        f"-fd {alpha} "
        f"-r {n_restarts} "
        "-e 0.20 "
        "-a "
        f"-l {n_iters} "
        f"-o {tmpdir.name}/scite.output > {tmpdir.name}/scite.log"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    with open(f"{tmpdir.name}/scite.output_ml0.gv") as fin:
        with open(f"{tmpdir.name}/scite.output_ml0_quoted.gv", "w") as fout:
            for line in fin:
                if " -> " in line:
                    line = line.strip()
                    a = line.split(" -> ")[0]
                    b = line.split(" -> ")[1].replace(";", "")
                    fout.write(f'"{a}" -> "{b}";\n')
                else:
                    fout.write(line)

    detail = {}
    with open(f"{tmpdir.name}/scite.log") as fin:
        for line in fin:
            line = line.strip()
            if "best value for beta:" in line:
                detail["beta"] = float(line.replace("best value for beta:", "").strip())
            if "best log score for tree:" in line:
                detail["score"] = float(
                    line.replace("best log score for tree:", "").strip()
                )

    G = nx.drawing.nx_pydot.read_dot(f"{tmpdir.name}/scite.output_ml0_quoted.gv")
    df_output = df_input.copy()
    for col in df_output.columns:
        df_output[col].values[:] = 0
    for i in range(df_output.shape[0]):
        muts = nx.shortest_path(G, source="Root", target=f"s{i}")
        muts.remove("Root")
        muts.remove(f"s{i}")
        if len(muts) > 0:
            df_output.loc[df_output.index[i], muts] = 1

    tmpdir.cleanup()

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
    experiment=False,
):
    tsc.logg.info(
        f"running infSCITE with alpha={alpha}, beta={beta}, n_iters={n_iters}, "
        f"n_restarts={n_restarts}"
    )
    tmpdir = tsc.ul.tmpdirsys(suffix=".infscite")

    np.savetxt(
        f"{tmpdir.name}/infscite.SC.T", df_input.values.T, delimiter="\t", fmt="%1.0f"
    )
    with open(f"{tmpdir.name}/infscite.geneNames", "w") as fout:
        fout.write("\n".join(df_input.columns))

    infscite = tsc.ul.get_file("trisicell.external/bin/infSCITE")
    cmd = (
        f"{infscite} "
        f"-i {tmpdir.name}/infscite.SC.T "
        f"-names {tmpdir.name}/infscite.geneNames "
        f"-n {df_input.shape[1]} "
        f"-m {df_input.shape[0]} "
        f"-ad {beta} "
        f"-fd {alpha} "
        f"-r {n_restarts} "
        "-z "
        "-e 0.20 "
        "-a "
        f"-l {n_iters} "
        f"-o {tmpdir.name}/infscite > {tmpdir.name}/infscite.log"
    )
    # "-rec 1 "
    # "-d "
    # "-s -p 10000 "

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    with open(f"{tmpdir.name}/infscite_ml0.gv") as fin:
        with open(f"{tmpdir.name}/infscite_ml0_quoted.gv", "w") as fout:
            for line in fin:
                if " -> " in line:
                    line = line.strip()
                    a = line.split(" -> ")[0]
                    b = line.split(" -> ")[1].replace(";", "")
                    fout.write(f'"{a}" -> "{b}";\n')
                else:
                    fout.write(line)

    detail = {}
    with open(f"{tmpdir.name}/infscite.log") as fin:
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

    G = nx.drawing.nx_pydot.read_dot(f"{tmpdir.name}/infscite_ml0_quoted.gv")
    df_output = df_input.copy()
    for col in df_output.columns:
        df_output[col].values[:] = 0
    for i in range(df_output.shape[0]):
        muts = nx.shortest_path(G, source=f"{df_input.shape[0]}", target=f"{i+1}")
        # muts.remove("Root")
        muts.remove(f"{i+1}")
        if len(muts) > 0:
            df_output.loc[df_output.index[i], muts] = 1

    tmpdir.cleanup()

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        for k, v in detail.items():
            tsc.logg.info(f"{k}: {v}")
        return df_output
    else:
        return df_output, running_time, detail["score"], detail["beta"]
