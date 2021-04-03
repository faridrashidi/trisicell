import itertools

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc


def combine_replicates_in_cfmatrix(df):
    df2 = df.groupby(df.index.str.split("_").str[0]).transform("prod")
    df2 = df2.groupby(df2.index.str.split("_").str[0]).first()
    return df2


def consensus_of_cfmatrices(df1, df2):
    df1_temp = df1.loc[df1.index]
    df2_temp = df2.loc[df1.index]
    mut1 = pd.DataFrame(df1.columns, columns=["muts"])
    mut1[["ensemble", "gene", "chrom", "position", "reference", "alteration"]] = (
        mut1["muts"].apply(tsc.ul.split_mut).apply(pd.Series)
    )
    mut2 = pd.DataFrame(df2.columns, columns=["muts"])
    mut2[["ensemble", "gene", "chrom", "position", "reference", "alteration"]] = (
        mut2["muts"].apply(tsc.ul.split_mut).apply(pd.Series)
    )
    mut1 = mut1.set_index("muts")
    mut2 = mut2.set_index("muts")
    mut1_t = (
        mut1["chrom"].astype(str)
        + ":"
        + mut1["position"].astype(str)
        + ":"
        + mut1["reference"].astype(str)
        + ":"
        + mut1["alteration"].astype(str)
    )
    mut2_t = (
        mut2["chrom"].astype(str)
        + ":"
        + mut2["position"].astype(str)
        + ":"
        + mut2["reference"].astype(str)
        + ":"
        + mut2["alteration"].astype(str)
    )
    mut1_t = pd.Series(mut1_t.index, index=mut1_t)
    mut2_t = pd.Series(mut2_t.index, index=mut2_t)
    shared = np.intersect1d(mut1_t.index, mut2_t.index)
    finals = []
    for mut1, mut2 in zip(mut1_t[shared].values, mut2_t[shared].values):
        if df1_temp[mut1].equals(df2_temp[mut2]):
            finals.append(mut1)
    return df1[finals]


def consensus_run(sc1, sc2):
    """Building the consensus tree between to phylogenetic trees.

    Parameters
    ----------
    df1 : :class:`pandas.DataFrame`
        First conflict-free matrix.
    df2 : :class:`pandas.DataFrame`
        Second conflict-free matrix.

    Returns
    -------
    :class:`networkx.DiGraph`
        Two trees derived by building the consensus procedure.
    """

    def get_cnt_tree(tree):
        cnt_tree = tree.copy()
        for u, v, l in cnt_tree.edges.data("label"):
            l = len(l.split(cnt_tree.graph["splitter_mut"]))
            cnt_tree.add_edge(u, v, label=l)
            if cnt_tree.in_degree(v) != 0:
                if "––" not in tree.nodes[v]["label"]:
                    cnt_tree.nodes[v]["label"] = cnt_tree.nodes[v]["label"].split(
                        cnt_tree.graph["splitter_cell"]
                    )
        return cnt_tree

    def get_node_labels(tree):
        nl = {}
        for x in tree.nodes:
            if tree.in_degree(x) != 0:
                if "––" not in tree.nodes[x]["label"]:
                    for y in tree.nodes[x]["label"]:
                        nl[y] = x
            else:
                nl["root"] = x
        return nl

    def has_path(tree, nodes, l1, l2):
        try:
            n1 = nodes[l1]
            n2 = nodes[l2]
            dij = nx.bidirectional_dijkstra(tree, n1, n2, weight="label")
            return True, dij[0]
        except:
            return False, None
        return False, None

    def merge_with_parents(tree, nodes, n1, n2):
        lca = nx.lowest_common_ancestor(tree, n1, n2)
        costs, paths = nx.bidirectional_dijkstra(tree, lca, n1, weight="label")
        paths = paths[::-1]
        labels = []
        for p in paths[1:]:
            # merge-with-parent, n1 -> n2 (n2 is going to be replaced by n1)
            if "––" not in tree.nodes[p]["label"]:
                if "––" not in tree.nodes[n1]["label"]:
                    tree.nodes[n1]["label"] += tree.nodes[p]["label"]
                else:
                    tree.nodes[n1]["label"] = tree.nodes[p]["label"]
                for x in tree.nodes[p]["label"]:
                    nodes[x] = n1
            cost = tree.edges[p, n1]["label"]
            tree = nx.contracted_nodes(tree, n1, p, self_loops=False)
        return tree, nodes, costs

    def _get_labels(cnt_tree):
        labels = {}
        [y for y in cnt_tree.nodes if "––" in cnt_tree.nodes[y]["label"]]
        for x in cnt_tree.nodes:
            if "root" not in cnt_tree.nodes[x]["label"]:
                lbs = []
                sub = nx.subgraph(
                    cnt_tree,
                    nx.algorithms.traversal.depth_first_search.dfs_tree(
                        cnt_tree, x
                    ).nodes,
                )
                for y in sub.nodes:
                    if "––" not in cnt_tree.nodes[y]["label"]:
                        lbs += cnt_tree.nodes[y]["label"]
                labels[x] = lbs
        return labels

    sc1_c = sc1.loc[:, sc1.sum() > 1].copy()
    sc2_c = sc2.loc[:, sc2.sum() > 1].copy()

    tree1 = tsc.ul.to_tree(sc1_c)
    tree2 = tsc.ul.to_tree(sc2_c)

    cnt_tree1 = get_cnt_tree(tree1)
    cnt_tree2 = get_cnt_tree(tree2)

    nodes1 = get_node_labels(cnt_tree1)
    nodes2 = get_node_labels(cnt_tree2)

    common_cells = list(np.intersect1d(nodes1.keys(), nodes2.keys())[0])

    total_cost = 0

    #### step 1,2
    for x, y in itertools.permutations(common_cells, 2):
        hp1, cost1 = has_path(cnt_tree1, nodes1, x, y)
        if hp1:
            hp2, cost2 = has_path(cnt_tree2, nodes2, x, y)
            if not hp2:
                tsc.logg.info("tree2:", x, y, cost2)
                cnt_tree2, nodes2, cost2 = merge_with_parents(
                    cnt_tree2, nodes2, nodes2[x], nodes2[y]
                )
                total_cost += cost2

    for x, y in itertools.permutations(common_cells, 2):
        hp2, cost2 = has_path(cnt_tree2, nodes2, x, y)
        if hp2:
            hp1, cost1 = has_path(cnt_tree1, nodes1, x, y)
            if not hp1:
                tsc.logg.info("tree1:", x, y, cost1)
                cnt_tree1, nodes1, cost1 = merge_with_parents(
                    cnt_tree1, nodes1, nodes1[x], nodes1[y]
                )
                total_cost += cost1

    for x, y in itertools.permutations(common_cells, 2):
        hp1, cost1 = has_path(cnt_tree1, nodes1, x, y)
        if hp1:
            hp2, cost2 = has_path(cnt_tree2, nodes2, x, y)
            if not hp2:
                tsc.logg.info("tree2:", x, y, cost2)
                cnt_tree2, nodes2, cost2 = merge_with_parents(
                    cnt_tree2, nodes2, nodes2[x], nodes2[y]
                )
                total_cost += cost2

    #### step 3
    labels1 = _get_labels(cnt_tree1)
    labels2 = _get_labels(cnt_tree2)

    matched = []
    for x, y in labels1.items():
        for p, q in labels2.items():
            if set(y) == set(q):
                matched.append((x, p))
    for i, j in matched:
        del labels1[i]
        del labels2[j]

    for y, _ in labels1.items():
        parent = list(cnt_tree1.predecessors(y))[0]
        tsc.logg.info("tree1:", "internal", y)
        cnt_tree1, nodes1, cost1 = merge_with_parents(cnt_tree1, nodes1, y, parent)
        total_cost += cost1

    for y, _ in labels2.items():
        parent = list(cnt_tree2.predecessors(y))[0]
        tsc.logg.info("tree2:", "internal", y)
        cnt_tree2, nodes2, cost2 = merge_with_parents(cnt_tree2, nodes2, y, parent)
        total_cost += cost2

    tsc.logg.info("total:", total_cost)

    return cnt_tree1, cnt_tree2
