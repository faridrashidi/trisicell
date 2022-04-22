import itertools

import ete3
import networkx as nx
import numpy as np

import trisicell as tsc
from trisicell.external._mp3 import build_tree, similarity
from trisicell.ul._trees import _to_newick


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


def PCSS(df_grnd, df_sol):
    """Pairwise cell shortest-path similarity score.

    For every pair of cells :math:`i` and :math:`j`, we computed the shortest-path
    :math:`d_{ij}` between the two cells in each tree. If the two cells belong to the
    same clone, their shortest-path distance is 0, otherwise the shortest-path distance
    equals the number of edges (regardless of direction) that separate the clones of the
    two cells. Finally, we summed up the absolute differences between the shortest-path
    distances of all unordered pairs of cells in the two trees.
    This measure is metric. The proof is given in the paper.

    This measure was introduced in :cite:`OncoNEM`.

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


def rf(df_grnd, df_sol):
    """Robinson-Foulds score.

    The Robinson–Foulds or symmetric difference metric is defined as (A + B) where A is
    the number of partitions of data implied by the first tree but not the second tree
    and B is the number of partitions of data implied by the second tree but not the
    first tree (although some software implementations divide the RF metric by 2
    and others scale the RF distance to have a maximum value of 1).

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

    inter = np.intersect1d(df_grnd.index, df_sol.index)
    if len(inter) == 0:
        tsc.logg.error("No common cells found between two trees!")
    df_grnd1 = df_grnd.loc[inter]
    df_sol1 = df_sol.loc[inter]

    tree_grnd = tsc.ul.to_tree(df_grnd1)
    tree_sol = tsc.ul.to_tree(df_sol1)

    nwk_grnd = _to_newick(tree_grnd)
    nwk_sol = _to_newick(tree_sol)

    # from skbio import TreeNode
    # tree_grnd = TreeNode.read([nwk_grnd])
    # tree_sol = TreeNode.read([nwk_sol])
    # return tree_grnd.compare_rfd(tree_sol)

    tree_grnd = ete3.Tree(nwk_grnd, format=1)
    tree_sol = ete3.Tree(nwk_sol, format=1)

    rf = tree_grnd.robinson_foulds(tree_sol, unrooted_trees=True)
    return 1 - rf[0] / rf[1]
