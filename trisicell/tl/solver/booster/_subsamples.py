import contextlib
import functools
import multiprocessing
import os
import time

import joblib
import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph
from tqdm import tqdm

import trisicell as tsc


def with_timeout(timeout):
    def decorator(decorated):
        @functools.wraps(decorated)
        def inner(*args, **kwargs):
            pool = multiprocessing.pool.ThreadPool(1)
            async_result = pool.apply_async(decorated, args, kwargs)
            try:
                return async_result.get(timeout)
            except multiprocessing.TimeoutError:
                return None

        return inner

    return decorator


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


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
    tmpdir,
    disable_tqdm,
):

    # tsc.logg.info(
    #     f"running Booster with alpha={alpha}, beta={beta}, solver={solver},"
    #     f" sample_on={sample_on}, sample_size={sample_size}, n_samples={n_samples},"
    #     f" n_jobs={n_jobs}, time_out={time_out}"
    # )
    # tsc.logg.info(f"id,m,i0,i1,i3,o0,o1,f01,f10,f30,f31,r")

    # graph = nx.DiGraph()
    # graph.add_nodes_from(df_input.columns)

    # def add_edge(muty, mutx):
    #     if graph.has_edge(muty, mutx):
    #         graph[muty][mutx]["label"] += 1
    #     else:
    #         graph.add_edge(muty, mutx, label=1)

    # def remove_edge(muty, mutx):
    #     if graph.has_edge(muty, mutx):
    #         graph[muty][mutx]["label"] -= 1

    @with_timeout(time_out)
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
                dfn, alpha, beta, time_out, 1, False, experiment=True
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

    s_time = time.time()
    with tqdm_joblib(
        tqdm(
            ascii=True,
            ncols=100,
            desc="SUBSAMPLING ",
            total=n_samples,
            position=0,
            disable=disable_tqdm,
        )
    ) as progress_bar:
        output = Parallel(n_jobs=n_jobs)(
            delayed(run)(i) for i in range(begin_sample, begin_sample + n_samples)
        )
    e_time = time.time()
    running_time = e_time - s_time

    # muts = np.append(["ROOT"], df_input.columns)
    # A = pd.DataFrame(0, index=muts, columns=muts)
    # D = pd.DataFrame(0, index=muts, columns=muts)
    # for i in range(n_samples):
    #     try:
    #         dfo = pd.read_csv(f"{tmpdir}/{i}.CFMatrix", sep="\t", index_col=0)
    #     except:
    #         continue
    #     for i in range(dfo.shape[1]):
    #         for j in range(dfo.shape[1]):
    #             if i != j:
    #                 mutx = dfo.columns[i]
    #                 muty = dfo.columns[j]
    #                 if (dfo[mutx] == dfo[muty]).all():
    #                     add_edge(muty, mutx)
    #                     A.loc[muty, mutx] += 1
    #                 elif (dfo[mutx] <= dfo[muty]).all():
    #                     add_edge(muty, mutx)
    #                     A.loc[muty, mutx] += 1
    #                 else:
    #                     D.loc[mutx, muty] += 1

    # for i in range(n_samples):
    #     try:
    #         dfo = pd.read_csv(f"{tmpdir}/{i}.CFMatrix", sep="\t", index_col=0)
    #     except:
    #         continue
    #     for i in range(dfo.shape[1]):
    #         for j in range(dfo.shape[1]):
    #             if i != j:
    #                 mutx = dfo.columns[i]
    #                 muty = dfo.columns[j]
    #                 if (dfo[mutx] == dfo[muty]).all():
    #                     pass
    #                 elif (dfo[mutx] <= dfo[muty]).all():
    #                     pass
    #                 else:
    #                     remove_edge(muty, mutx)
    #                     remove_edge(mutx, muty)

    # A.loc["ROOT", :] = np.inf
    # D.loc["ROOT", :] = 0
    # D.loc[:, "ROOT"] = 0
    # A.loc[:, "ROOT"] = 0

    # graph.add_nodes_from(["ROOT"])
    # for mut in df_input.columns:
    #     add_edge("ROOT", mut)

    # T = nx.algorithms.tree.maximum_spanning_arborescence(graph, attr="label")
    # mygraph = to_agraph(T)
    # mygraph.layout(prog="dot")
    # mygraph.draw(f"{tmpdir}/_tree.png")

    # tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
