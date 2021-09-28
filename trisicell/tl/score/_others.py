import itertools

import networkx as nx
import numpy as np

import trisicell as tsc
from trisicell.tl.score._mp3 import build_tree, similarity


def bourque(df_grnd, df_sol):
    """Bourque distances for mutation trees.

    This measure was introduced in :cite:`Bourque`.

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
        Similarity out of one.
    """

    # TODO: implement
    return None


def mp3(df_grnd, df_sol):
    """Triplet-based similarity score.

    For fully multilabeled trees with poly-occurring labels.
    This measure was introduced in :cite:`MP3`.

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
        Similarity out of one.
    """

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    if len(inter) == 0:
        tsc.logg.error("No common mutations found between two trees!")
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
    tree_grnd = tsc.ul.to_mtree(tree_grnd)
    tree_sol = tsc.ul.to_mtree(tree_sol)

    for n in tree_grnd.nodes:
        if tree_grnd.in_degree(n) > 0:
            tree_grnd.nodes[n]["label"] = ",".join(tree_grnd.nodes[n]["label"])
    for n in tree_sol.nodes:
        if tree_sol.in_degree(n) > 0:
            tree_sol.nodes[n]["label"] = ",".join(tree_sol.nodes[n]["label"])

    T1 = build_tree(tree_grnd)
    T2 = build_tree(tree_sol)
    return similarity(T1, T2)


def caset(df_grnd, df_sol):
    """Commonly Ancestor Sets score.

    This measure was introduced in :cite:`CASet_DISC`.

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
        Similarity out of one.
    """

    def _get_ancesteral_set(tree):
        root = tsc.ul.root_id(tree)
        ancesteral_set = {}
        for n in tree.nodes:
            if tree.in_degree(n) > 0 and "––" not in tree.nodes[n]["label"]:
                for mut in tree.nodes[n]["label"]:
                    ancester_set = ["root"]
                    ancester_set += tree.nodes[n]["label"]  # self ancestor is ok
                    for m in nx.shortest_path(tree, root, n):
                        if tree.in_degree(m) > 0 and "––" not in tree.nodes[m]["label"]:
                            ancester_set += tree.nodes[m]["label"]
                    ancesteral_set[mut] = ancester_set
        return ancesteral_set

    def _get_common_ancesteral_set(muts, ancesteral_set):
        common_ancesteral_set = {}
        for x, y in itertools.combinations(muts, 2):
            common_ancesteral_set[(x, y)] = np.intersect1d(
                ancesteral_set[x], ancesteral_set[y]
            )
        return common_ancesteral_set

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    if len(inter) == 0:
        tsc.logg.error("No common mutations found between two trees!")
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
    tree_grnd = tsc.ul.to_mtree(tree_grnd)
    tree_sol = tsc.ul.to_mtree(tree_sol)

    ancesteral_set_grnd = _get_ancesteral_set(tree_grnd)
    ancesteral_set_sol = _get_ancesteral_set(tree_sol)

    common_ancesteral_set_grnd = _get_common_ancesteral_set(inter, ancesteral_set_grnd)
    common_ancesteral_set_sol = _get_common_ancesteral_set(inter, ancesteral_set_sol)

    final = []
    for x, y in itertools.combinations(inter, 2):
        a = len(
            np.intersect1d(
                common_ancesteral_set_grnd[(x, y)], common_ancesteral_set_sol[(x, y)]
            )
        )
        b = len(
            np.union1d(
                common_ancesteral_set_grnd[(x, y)], common_ancesteral_set_sol[(x, y)]
            )
        )
        final.append(a / b)

    return np.mean(final)


def disc(df_grnd, df_sol):
    """Distinctly Inherited Sets score.

    This measure was introduced in :cite:`CASet_DISC`.

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
        Similarity out of one.
    """

    def _get_ancesteral_set(tree):
        root = tsc.ul.root_id(tree)
        ancesteral_set = {}
        for n in tree.nodes:
            if tree.in_degree(n) > 0 and "––" not in tree.nodes[n]["label"]:
                for mut in tree.nodes[n]["label"]:
                    ancester_set = ["root"]
                    ancester_set += tree.nodes[n]["label"]  # self ancestor is ok
                    for m in nx.shortest_path(tree, root, n):
                        if tree.in_degree(m) > 0 and "––" not in tree.nodes[m]["label"]:
                            ancester_set += tree.nodes[m]["label"]
                    ancesteral_set[mut] = ancester_set
        return ancesteral_set

    def _get_distinctly_inherited_set(muts, ancesteral_set):
        common_ancesteral_set = {}
        for x, y in itertools.permutations(muts, 2):
            common_ancesteral_set[(x, y)] = np.setdiff1d(
                ancesteral_set[x], ancesteral_set[y]
            )
        return common_ancesteral_set

    inter = np.intersect1d(df_grnd.columns, df_sol.columns)
    if len(inter) == 0:
        tsc.logg.error("No common mutations found between two trees!")
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
    tree_grnd = tsc.ul.to_mtree(tree_grnd)
    tree_sol = tsc.ul.to_mtree(tree_sol)

    ancesteral_set_grnd = _get_ancesteral_set(tree_grnd)
    ancesteral_set_sol = _get_ancesteral_set(tree_sol)

    distinctly_inherited_set_grnd = _get_distinctly_inherited_set(
        inter, ancesteral_set_grnd
    )
    distinctly_inherited_set_sol = _get_distinctly_inherited_set(
        inter, ancesteral_set_sol
    )

    final = []
    for x, y in itertools.combinations(inter, 2):
        a = len(
            np.intersect1d(
                distinctly_inherited_set_grnd[(x, y)],
                distinctly_inherited_set_sol[(x, y)],
            )
        )
        b = len(
            np.union1d(
                distinctly_inherited_set_grnd[(x, y)],
                distinctly_inherited_set_sol[(x, y)],
            )
        )
        if (
            b > 0
        ):  # FIXME: if a and b are in the same node distinctly_inherited_set is empty
            final.append(a / b)

    return np.mean(final)
