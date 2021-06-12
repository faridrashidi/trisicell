import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

import trisicell as tsc


def subsampling(
    df_input,
    alpha,
    beta,
    solver,
    sample_on,
    sample_size,
    n_samples,
    begin_sample,
    n_jobs,
    time_out,
    n_iterations,
    tmpdir,
    disable_tqdm,
):

    # tsc.logg.info(
    #     f"running Booster with alpha={alpha}, beta={beta}, solver={solver},"
    #     f" sample_on={sample_on}, sample_size={sample_size}, n_samples={n_samples},"
    #     f" n_jobs={n_jobs}, time_out={time_out}"
    # )
    # tsc.logg.info(f"id,m,i0,i1,i3,o0,o1,f01,f10,f30,f31,r")

    @tsc.ul.with_timeout(time_out)
    def run(i):
        if sample_on == "muts":
            dfn = df_input.sample(n=sample_size, replace=False, axis=1)
        elif sample_on == "cells":
            dfn = df_input.sample(n=sample_size, replace=False, axis=0)

        if solver.lower() == "phiscs":
            dfo = tsc.tl.phiscsb(dfn, alpha, beta, experiment=True)
            dfo.to_csv(f"{tmpdir}/{i}.CFMatrix", sep="\t")
        elif solver.lower() == "scite":
            dfo, _, _, _ = tsc.tl.scite(
                dfn, alpha, beta, n_iters=n_iterations, n_restarts=1, experiment=True
            )
            dfo.to_csv(f"{tmpdir}/{i}.CFMatrix", sep="\t")
        elif solver.lower() == "scistree":
            x = ((dfn == 1).sum() >= 2) & ((dfn != 3).sum() >= 0.2 * dfn.shape[0])
            dfn = dfn[dfn.columns[x]]
            if dfn.shape[1] < 2:
                return None
            dfo, run_time = tsc.tl.scistree(dfn, alpha, beta, False, experiment=True)
            i0 = np.sum(dfn.values == 0)
            i1 = np.sum(dfn.values == 1)
            i3 = np.sum(dfn.values == 3)
            o0 = np.sum(dfo.values == 0)
            o1 = np.sum(dfo.values == 1)
            f01, f10, f30, f31 = tsc.ul.count_flips(dfn.values, dfo.values)
            # tsc.logg.info(
            #     f"{i},{dfn.shape[1]},{i0},{i1},{i3},{o0},{o1},{f01},{f10},{f30},{f31},{run_time:.1f}"
            # )
            dfo.to_csv(f"{tmpdir}/{i}.CFMatrix", sep="\t")

    with tsc.ul.tqdm_joblib(
        tqdm(
            ascii=True,
            ncols=100,
            desc="SUBSAMPLING    (1/3)",
            total=n_samples,
            position=0,
            disable=disable_tqdm,
        )
    ) as progress_bar:
        output = Parallel(n_jobs=n_jobs)(
            delayed(run)(i) for i in range(begin_sample, begin_sample + n_samples)
        )
