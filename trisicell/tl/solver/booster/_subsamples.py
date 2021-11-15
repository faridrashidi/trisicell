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
    time_limit,
    n_iterations,
    tmpdir,
    disable_tqdm,
):

    if solver.lower() == "scite":
        time_limit = n_iterations

    @tsc.ul.with_timeout(time_limit)
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
            dfo, _ = tsc.tl.scistree(dfn, alpha, beta, experiment=True)
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
    ):
        Parallel(n_jobs=n_jobs)(
            delayed(run)(i) for i in range(begin_sample, begin_sample + n_samples)
        )
