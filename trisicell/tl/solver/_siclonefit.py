import os
import time

import numpy as np
import pandas as pd

import trisicell as tsc


def siclonefit(df_input, alpha, beta, n_iters):
    # TODO: implement
    executable = tsc.ul.executable("SiCloneFiTComplete.jar", "SiCloneFit")

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
    I_mtr = df_input.values

    cmd = (
        f"java -jar {executable} "
        f"-m {df_input.shape[0]} "
        f"-n {df_input.shape[1]} "
        f"-ipMat {tmpdir.name}/siclonefit.input "
        f"-fp {alpha} "
        f"-fn {beta} "
        # "-df 0 "
        f"-missing {np.sum(I_mtr == 3)/(I_mtr.size)} "
        # "-f 3 "
        # "-recurProb 0 "
        # "-delProb 0 "
        # "-LOHProb 0 "
        # f"-iter {n_iters} "
        f"-cellNames {tmpdir.name}/siclonefit.cellnames "
        f"-geneNames {tmpdir.name}/siclonefit.genenames "
        # "-r "
        # "-burnin "
        # "-printIter "
        # "-treeIter "
        # "-doublet "
        f"-outDir {tmpdir.name} > {tmpdir.name}/siclonefit.log"
    )
    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    df_output = pd.read_csv(
        f"{tmpdir.name}/samples/best/best_MAP_predicted_genotype.txt",
        sep=" ",
        header=None,
        index_col=0,
    ).T
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    tmpdir.cleanup()

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
