import networkx as nx
import numpy as np

import trisicell as tsc


def _get_cnt_tree(tree):
    cnt_tree = tree.copy()
    for u, v, _ in cnt_tree.edges.data("label"):
        cnt_tree.add_edge(u, v, label="")
        if cnt_tree.in_degree(v) != 0:
            if "––" not in tree.nodes[v]["label"]:
                cnt_tree.nodes[v]["label"] = cnt_tree.nodes[v]["label"].split(
                    cnt_tree.graph["splitter_cell"]
                )
    return cnt_tree


def _get_edge_set(tree):
    edges = {}
    for e in tree.edges():
        nodes = nx.algorithms.traversal.depth_first_search.dfs_tree(tree, e[1]).nodes
        cells = []
        for n in nodes:
            if "––" not in tree.nodes[n]["label"]:
                cells += tree.nodes[n]["label"]
        edges[e] = cells
    return edges


def consensus(sc1, sc2):
    """Build the consensus tree between two tumor progression trees.

    Parameters
    ----------
    df1 : :class:`pandas.DataFrame`
        First conflict-free matrix.
    df2 : :class:`pandas.DataFrame`
        Second conflict-free matrix.

    Returns
    -------
    :class:`networkx.DiGraph`
        The consensus tree in which cells are in edges and mutations are in edges.
    """

    good_cells1 = sc1.index[(sc1 != 0).sum(axis=1) > 0]
    good_cells2 = sc2.index[(sc2 != 0).sum(axis=1) > 0]
    common_cells = np.intersect1d(good_cells1, good_cells2)
    if len(common_cells) == 0:
        tsc.logg.error("No common cells found in two inputs!")
    sc1_p = sc1.loc[common_cells, :]
    sc2_p = sc2.loc[common_cells, :]

    tree1 = tsc.ul.to_tree(sc1_p)
    tree2 = tsc.ul.to_tree(sc2_p)

    cnt_tree1 = _get_cnt_tree(tree1)
    cnt_tree2 = _get_cnt_tree(tree2)

    edge_set1 = _get_edge_set(cnt_tree1)
    edge_set2 = _get_edge_set(cnt_tree2)

    return cnt_tree1, cnt_tree2, edge_set1, edge_set2
