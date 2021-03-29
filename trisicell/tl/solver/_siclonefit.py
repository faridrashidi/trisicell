import os

import numpy as np
import pandas as pd

import trisicell as tsc


def siclonefit(df_input, alpha, beta, n_iters, save_inter=False):
    tsc.logg.info(
        f"running SiCloneFit with alpha={alpha}, beta={beta}, n_iters={n_iters}"
    )
    tmpdir = tsc.ul.tmpdir(prefix="trisicell.", suffix=".siclonefit", dirname=".")

    df_input.T.reset_index(drop=True).to_csv(
        f"{tmpdir}/siclonefit.input", sep=" ", header=None
    )
    with open(f"{tmpdir}/siclonefit.cellnames", "w") as fout:
        fout.write(" ".join(df_input.index))
    with open(f"{tmpdir}/siclonefit.genenames", "w") as fout:
        fout.write(" ".join(df_input.columns))
    I = df_input.values

    siclonefit = tsc.ul.get_file("trisicell.external/bin/SiCloneFiTComplete.jar")
    cmd = (
        f"java -jar {siclonefit} "
        f"-m {df_input.shape[0]} "
        f"-n {df_input.shape[1]} "
        f"-ipMat {tmpdir}/siclonefit.input "
        f"-fp {alpha} "
        f"-fn {beta} "
        "-df 0 "
        f"-missing {np.sum(I == 3)/(I.size)} "
        "-f 3 "
        "-recurProb 0 "
        "-delProb 0 "
        "-LOHProb 0 "
        f"-iter {n_iters} "
        f"-cellNames {tmpdir}/siclonefit.cellnames "
        f"-geneNames {tmpdir}/siclonefit.genenames "
        f"-outDir {tmpdir} > {tmpdir}/siclonefit.log"
    )
    # check the following parameters
    # -burnin
    # -printIter
    # -treeIter
    # -doublet
    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    df = pd.read_csv(
        f"{tmpdir}/20p_missing_samples/best/best_MAP_predicted_genotype.txt",
        sep=" ",
        header=None,
        index_col=0,
    ).T
    df_output = pd.DataFrame(df.values)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    if not save_inter:
        tsc.ul.cleanup(tmpdir)

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
