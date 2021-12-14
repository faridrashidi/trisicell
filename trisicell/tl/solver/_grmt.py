import os
import time

import networkx as nx

import trisicell as tsc


def grmt(df_input, alpha, beta, n_iters=30, n_threads=1):
    """Solving using GRMT.

    Generative Reconstruction of Mutation Tree From Scratch Using Single-Cell
    Sequencing Data :cite:`GRMT`.

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
    n_iters : :obj:`int`, optional
        Number of iterations, by default 30
    n_threads : :obj:`int`, optional
        Number of threads, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    executable = tsc.ul.executable("grmt", "GRMT")

    tsc.logg.info(
        f"running GRMT with alpha={alpha}, beta={beta}, n_threads={n_threads}"
    )

    tmpdir = tsc.ul.tmpdirsys(suffix=".sphyr")

    df_input.to_csv(f"{tmpdir.name}/grmt.input", sep="\t", header=None, index=None)
    with open(f"{tmpdir.name}/grmt.cellnames", "w") as fout:
        fout.write("\n".join(df_input.index) + "\n")
    with open(f"{tmpdir.name}/grmt.mutnames", "w") as fout:
        fout.write("\n".join(df_input.columns) + "\n")

    cmd = (
        f"{executable} "
        f"--input {tmpdir.name}/grmt.input "
        f"--alpha {alpha} "
        f"--beta {beta} "
        f"--clabel {tmpdir.name}/grmt.cellnames "
        f"--mlabel {tmpdir.name}/grmt.mutnames "
        "--lambda 0.7 "
        "--kappa 1 "
        "--n_init 100 "
        f"--n_iter {n_iters} "
        "--maxl 0 "
        f"--threads {n_threads} "
        f"--output {tmpdir.name}/grmt "
        f"&> {tmpdir.name}/grmt.log"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    df_output = df_input.copy()
    df_output[:] = 0
    G = nx.drawing.nx_pydot.read_dot(f"{tmpdir.name}/grmt.dot")
    relation = {
        v.replace('"', "").replace("\\n", ""): k
        for k, v in nx.get_node_attributes(G, "label").items()
    }
    data = {}
    for k, v in relation.items():
        if " " in k:
            for t in k.split(" "):
                data[t] = v
        else:
            data[k] = v
    for cell in df_output.index:
        muts = nx.shortest_path(G, source="0", target=data[cell])
        muts.remove("0")
        muts.remove(data[cell])
        muts = [G.nodes[m]["label"].replace('"', "") for m in muts]
        if len(muts) > 0:
            df_output.loc[cell, muts] = 1

    tmpdir.cleanup()

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
