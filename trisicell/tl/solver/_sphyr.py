import os
import time

import pandas as pd

import trisicell as tsc


def sphyr(
    df_input,
    alpha,
    beta,
    n_restarts=10,
    n_threads=1,
    time_limit=None,
    n_cell_clusters=10,
    n_mut_clusters=15,
):
    """Solving using SPhyR.

    Tumor phylogeny estimation from single-cell sequencing data under loss and error
    :cite:`SPhyR`.

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
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    # TODO: implement
    if not os.path.exists(f"{tsc.settings.tools}/kDPFC"):
        tsc.logg.error("Cannot find the binary file of SPhyR with `kDPFC` name!")

    tsc.logg.info(
        f"running SPhyR with alpha={alpha}, beta={beta}, n_restarts={n_restarts}, "
        f"n_threads={n_threads}, time_limit={time_limit}, "
        f"n_cell_clusters={n_cell_clusters}, n_mut_clusters={n_mut_clusters}"
    )

    tmpdir = tsc.ul.tmpdirsys(suffix=".sphyr")
    # tmpdir = tsc.ul.tmpdir(suffix=".sphyr")
    with open(f"{tmpdir.name}/sphyr.input", "a") as fout:
        fout.write(f"{df_input.shape[0]} #cells\n{df_input.shape[1]} #SNVs\n")
        df_input.replace(3, -1).to_csv(fout, sep=" ", header=None, index=None)
    with open(f"{tmpdir.name}/sphyr.cellnames", "w") as fout:
        fout.write("\n".join(df_input.index) + "\n")
    with open(f"{tmpdir.name}/sphyr.genenames", "w") as fout:
        fout.write("\n".join(df_input.columns) + "\n")

    cmd = (
        f"{tsc.settings.tools}/kDPFC "
        f"{tmpdir.name}/sphyr.input "
        f"-a {alpha} "
        f"-b {beta} "
        f"-N {n_restarts} "
        f"-t {n_threads} "
        # f"-lC {n_mut_clusters} "
        # f"-lT {n_cell_clusters} "
        f"-T {time_limit if time_limit is not None else -1} "
        "-k 0 "
        f"> {tmpdir.name}/sphyr.output "
        f"2> {tmpdir.name}/sphyr.log"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    df_output = pd.read_csv(
        f"{tmpdir.name}/sphyr.output",
        sep=" ",
        skiprows=[0, 1],
        header=None,
    )
    df_output.index = pd.read_csv(f"{tmpdir.name}/sphyr.cellnames", header=None)[0]
    df_output.columns = pd.read_csv(f"{tmpdir.name}/sphyr.genenames", header=None)[0]
    df_output.index.name = "cellIDxmutID"

    tmpdir.cleanup()

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    # cmd = (
    #     f"{sphyr}/visualize "
    #     f"{tmpdir}/sphyr.output "
    #     f"-c {tmpdir}/sphyr.genenames "
    #     f"-t {tmpdir}/sphyr.cellnames "
    #     f"> {tmpdir}/sphyr.dot"
    # )
    # os.system(cmd)
    # cmd = f"dot -Tpng {tmpdir}/sphyr.dot -o {tmpdir}/sphyr.png"
    # os.system(cmd)

    return df_output
