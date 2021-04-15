import os
import time

import apted
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import trisicell as tsc
from trisicell.external._mltd import run_mltd
from trisicell.ul._trees import _split_labels, _to_apted


def ad(df_grnd, df_sol):
    """Ancestor-descendent accuracy.

    For each pair of mutations in ground truth tree that are in
    ancestor-descendant relation (same nodes excluded) we check
    whether this relationship is preserved in the inferred tree.

    This measure was introduced in :cite:`B-SCITE`.

    Parameters
    ----------
    df_grnd : :class:`pandas.DataFrame`
        The first genotype matrix (e.g. ground truth)
        This matrix must be conflict-free.
    df_sol : :class:`pandas.DataFrame`
        The second genotype matrix (e.g. solution/inferred)
        This matrix must be conflict-free.

    Returns
    -------
    :obj:`float`
        The score is out of 1.00
    """

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    M_grnd = df_grnd[inter].values
    M_sol = df_sol[inter].values
    error_pairs = []
    n_adpairs = 0
    for i in range(M_grnd.shape[1]):
        for j in range(i, M_grnd.shape[1]):
            cap1 = M_grnd[:, i] * M_grnd[:, j]
            cap2 = M_sol[:, i] * M_sol[:, j]
            if np.sum(cap1) > 0 and np.sum(M_grnd[:, i]) != np.sum(M_grnd[:, j]):
                n_adpairs = n_adpairs + 1
                if np.sum(cap2) == 0:
                    error_pairs.append([i, j])
                else:
                    if np.sum(M_grnd[:, j]) > np.sum(M_grnd[:, i]) and np.sum(
                        M_sol[:, j]
                    ) <= np.sum(M_sol[:, i]):
                        error_pairs.append([i, j])
                    else:
                        if np.sum(M_grnd[:, i]) > np.sum(M_grnd[:, j]) and np.sum(
                            M_sol[:, i]
                        ) <= np.sum(M_sol[:, j]):
                            error_pairs.append([i, j])
    return 1 - len(error_pairs) / n_adpairs


def dl(df_grnd, df_sol):
    """Different-lineage accuracy.

    For each pair of mutations in ground truth tree that are
    in different-lineages relation we check whether the same relationship
    is preserved in the inferred tree.

    This measure was introduced in :cite:`B-SCITE`.

    Parameters
    ----------
    df_grnd : :class:`pandas.DataFrame`
        The first genotype matrix (e.g. ground truth)
        This matrix must be conflict-free.
    df_sol : :class:`pandas.DataFrame`
        The second genotype matrix (e.g. solution/inferred)
        This matrix must be conflict-free.

    Returns
    -------
    :obj:`float`
        The score is out of 1.00
    """

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    M_grnd = df_grnd[inter].values
    M_sol = df_sol[inter].values
    error_pairs = []
    n_dlpairs1 = 0
    n_dlpairs2 = 0
    for i in range(M_grnd.shape[1]):
        for j in range(i, M_grnd.shape[1]):
            cap1 = M_grnd[:, i] * M_grnd[:, j]
            cap2 = M_sol[:, i] * M_sol[:, j]
            if (
                np.sum(cap1) == 0
                and np.sum(M_grnd[:, i]) != 0
                and np.sum(M_grnd[:, j]) != 0
            ):
                n_dlpairs1 = n_dlpairs1 + 1
                if (
                    np.sum(cap2) == 0
                    and np.sum(M_sol[:, i]) != 0
                    and np.sum(M_sol[:, j]) != 0
                ):
                    n_dlpairs2 = n_dlpairs2 + 1
    return n_dlpairs2 / n_dlpairs1


def cc(df_grnd, df_sol):
    """Co-clustering accuracy.

    For each pair of mutations in ground truth tree that are on the same node we look
    relationship is preserved in the inferred tree.

    This measure was introduced in :cite:`B-SCITE`.

    Parameters
    ----------
    df_grnd : :class:`pandas.DataFrame`
        The first genotype matrix (e.g. ground truth)
        This matrix must be conflict-free.
    df_sol : :class:`pandas.DataFrame`
        The second genotype matrix (e.g. solution/inferred)
        This matrix must be conflict-free.

    Returns
    -------
    :obj:`float`
        The score is out of 1.00
    """

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    M_grnd = df_grnd[inter].values
    M_sol = df_sol[inter].values

    # TODO: implement.
    pass


def mltd(df_grnd, df_sol):
    """Multi-labeled tree dissimilarity measure (MLTD).

    This measure was introduced in :cite:`MLTD`.

    Parameters
    ----------
    df_grnd : :class:`pandas.DataFrame`
        The first genotype matrix (e.g. ground truth)
        This matrix must be conflict-free.
    df_sol : :class:`pandas.DataFrame`
        The second genotype matrix (e.g. solution/inferred)
        This matrix must be conflict-free.

    Returns
    -------
    :obj:`dict`
        {'distance', 'similarity', 'normalized_similarity'}
    """

    def _convert_tree_to_mtld_input(tree, file):
        with open(file, "w") as fout:
            for u, v, l in tree.edges.data("label"):
                if tree.in_degree(u) == 0:
                    fout.write(f"{u}=\n")
                muts = l.split(tree.graph["splitter_mut"])
                fout.write(f"{v}={','.join(muts)}\n")
            for u in tree.nodes:
                children = [str(n) for n in tree.neighbors(u)]
                if len(children) == 0:
                    continue
                fout.write(f"{u}:{','.join(children)}\n")

    tmpdir = tsc.ul.tmpdirsys(suffix=".mltd")

    df_grnd.columns = df_grnd.columns.str.replace(":", "_").str.replace("=", "_")
    df_sol.columns = df_sol.columns.str.replace(":", "_").str.replace("=", "_")
    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    df_grnd1 = df_grnd[inter]
    df_sol1 = df_sol[inter]

    tree_grnd = tsc.ul.to_tree(df_grnd1)
    tree_sol = tsc.ul.to_tree(df_sol1)

    _convert_tree_to_mtld_input(tree_grnd, f"{tmpdir.name}/grnd.in")
    _convert_tree_to_mtld_input(tree_sol, f"{tmpdir.name}/sol.in")

    s_time = time.time()
    result = run_mltd(f"{tmpdir.name}/grnd.in", f"{tmpdir.name}/sol.in")
    e_time = time.time()
    running_time = e_time - s_time

    tmpdir.cleanup()

    return result


def tpted(df_grnd, df_sol):
    """Tumor phylogeny tree edit distance measure (TPTED).

    This measure was introduced in :cite:`PhISCS`. This implementation uses
    `APTED <https://github.com/JoaoFelipe/apted>`_.

    Parameters
    ----------
    df_grnd : :class:`pandas.DataFrame`
        The first genotype matrix (e.g. ground truth)
        This matrix must be conflict-free.
    df_sol : :class:`pandas.DataFrame`
        The second genotype matrix (e.g. solution/inferred)
        This matrix must be conflict-free.

    Returns
    -------
    :obj:`int`
        distance
    """

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    df_grnd1 = df_grnd[inter]
    df_sol1 = df_sol[inter]

    tree_grnd = tsc.ul.to_tree(df_grnd1)
    tree_sol = tsc.ul.to_tree(df_sol1)
    inter = np.setdiff1d(
        inter,
        np.union1d(
            tree_grnd.graph["become_germline"], tree_sol.graph["become_germline"]
        ),
    )

    df_grnd1 = df_grnd[inter]
    df_sol1 = df_sol[inter]

    tree_grnd = tsc.ul.to_tree(df_grnd1)
    tree_sol = tsc.ul.to_tree(df_sol1)

    mt_grnd = tsc.ul.to_mtree(tree_grnd)
    mt_sol = tsc.ul.to_mtree(tree_sol)

    sl_grnd = _split_labels(mt_grnd, mt_guide=mt_sol)
    sl_sol = _split_labels(mt_sol, mt_guide=mt_sol)

    apted_grnd = _to_apted(sl_grnd)
    apted_sol = _to_apted(sl_sol)

    tree1 = apted.helpers.Tree.from_text(apted_grnd)
    tree2 = apted.helpers.Tree.from_text(apted_sol)

    ap = apted.APTED(tree1, tree2)
    ed = ap.compute_edit_distance()

    return 1 - ed / (2 * len(inter))
