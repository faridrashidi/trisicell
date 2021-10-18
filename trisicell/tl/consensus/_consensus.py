import networkx as nx
import numpy as np

import trisicell as tsc


def _get_cnt_tree(tree):
    cnt_tree = tree.copy()
    for u, v, muts in cnt_tree.edges.data("label"):
        mutations = muts.split(cnt_tree.graph["splitter_mut"])
        cnt_tree.add_edge(u, v, label=len(mutations), mutations=mutations)
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
        edges[tuple(sorted(cells))] = e
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

    lists = list(edge_set1.items())
    lists.sort(key=lambda x: len(x[0]))
    final_tree = cnt_tree1.copy()
    cost = 0
    for c, e in lists:
        if c not in edge_set2:
            if final_tree.nodes[e[0]]["label"] == "––":
                label = final_tree.nodes[e[1]]["label"]
            else:
                if final_tree.nodes[e[1]]["label"] == "––":
                    label = final_tree.nodes[e[0]]["label"]
                else:
                    label = (
                        final_tree.nodes[e[0]]["label"]
                        + final_tree.nodes[e[1]]["label"]
                    )
            cost += final_tree.edges[e]["label"]
            final_tree = nx.contracted_nodes(final_tree, e[0], e[1], self_loops=False)
            final_tree.nodes[e[0]]["label"] = label
    tsc.logg.info("    ---> total cost:", cost)

    # convert back to the trisicell tree format
    for v in final_tree.nodes:
        if final_tree.in_degree(v) != 0:
            if "––" not in final_tree.nodes[v]["label"]:
                final_tree.nodes[v]["label"] = final_tree.graph["splitter_cell"].join(
                    final_tree.nodes[v]["label"]
                )
    i = 0
    for e, u, _ in final_tree.edges.data("label"):
        final_tree.edges[(e, u)]["label"] = f"m{i}"
        i += 1
    data = tsc.ul.to_cfmatrix(final_tree)
    final_tree = tsc.ul.to_tree(data)

    return final_tree
