import os

import numpy as np
import pandas as pd

import trisicell as tsc


def siclonefit(df_input, alpha, beta, n_iters):
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
    I = df_input.values

    siclonefit = tsc.ul.get_file("trisicell.external/bin/SiCloneFiTComplete.jar")
    cmd = (
        f"java -jar {siclonefit} "
        f"-m {df_input.shape[0]} "
        f"-n {df_input.shape[1]} "
        f"-ipMat {tmpdir.name}/siclonefit.input "
        f"-fp {alpha} "
        f"-fn {beta} "
        "-df 0 "
        f"-missing {np.sum(I == 3)/(I.size)} "
        "-f 3 "
        "-recurProb 0 "
        "-delProb 0 "
        "-LOHProb 0 "
        f"-iter {n_iters} "
        f"-cellNames {tmpdir.name}/siclonefit.cellnames "
        f"-geneNames {tmpdir.name}/siclonefit.genenames "
        f"-outDir {tmpdir.name} > {tmpdir.name}/siclonefit.log"
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
        f"{tmpdir.name}/20p_missing_samples/best/best_MAP_predicted_genotype.txt",
        sep=" ",
        header=None,
        index_col=0,
    ).T
    df_output = pd.DataFrame(df.values)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    tmpdir.cleanup()

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
