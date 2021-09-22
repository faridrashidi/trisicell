import networkx as nx
import numpy as np
import pandas as pd

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


def _add_artificial_edge(tree):
    artificial_added = []
    cp_tree = tree.copy()
    nodes = list(cp_tree.nodes)
    max_node = max(nodes)
    for n in nodes:
        if cp_tree.in_degree(n) == 1 and cp_tree.out_degree(n) > 0:
            if "––" not in cp_tree.nodes[n]["label"]:
                for x in cp_tree.nodes[n]["label"]:
                    artificial_added.append(x)
                    max_node += 1
                    cp_tree.add_node(max_node, label=[x])
                    cp_tree.add_edge(n, max_node)
                cp_tree.nodes[n]["label"] = "––"
        elif cp_tree.in_degree(n) == 1 and cp_tree.out_degree(n) == 0:
            if len(cp_tree.nodes[n]["label"]) > 1:
                for x in cp_tree.nodes[n]["label"]:
                    max_node += 1
                    cp_tree.add_node(max_node, label=[x])
                    cp_tree.add_edge(n, max_node)
                cp_tree.nodes[n]["label"] = "––"
    return cp_tree, artificial_added


def _find_set_of_internal_nodes(tree):
    result = []
    for n in tree.nodes:
        if tree.in_degree(n) == 1 and tree.out_degree(n) > 0:
            tmp = []
            nodes = nx.algorithms.traversal.depth_first_search.dfs_tree(tree, n).nodes
            for m in nodes:
                if "––" not in tree.nodes[m]["label"]:
                    tmp += tree.nodes[m]["label"]
            result.append(sorted(tmp))
    tmp = []
    for m in tree.nodes:
        if "––" not in tree.nodes[m]["label"] and tree.in_degree(m) == 1:
            tmp += tree.nodes[m]["label"]
    result.append(sorted(tmp))
    return result


def consensus_day(sc1, sc2):
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

    cnt_tree1, artificial_added1 = _add_artificial_edge(cnt_tree1)
    cnt_tree2, artificial_added2 = _add_artificial_edge(cnt_tree2)
    artificial_added = list(np.union1d(artificial_added1, artificial_added2))
    tsc.logg.debug(artificial_added)

    soin1 = _find_set_of_internal_nodes(cnt_tree1)
    soin2 = _find_set_of_internal_nodes(cnt_tree2)

    intersection = list(np.intersect1d(soin1, soin2))
    intersection_len = list(map(len, intersection))
    valmax = max(intersection_len)
    argmax = intersection_len.index(valmax)
    indices = intersection[argmax]
    columns = [f"m{x}" for x in range(len(intersection))]
    df = pd.DataFrame(0, index=indices, columns=columns, dtype=int)
    for i, inter in enumerate(intersection):
        df.loc[inter, f"m{i}"] = 1

    final_tree = tsc.ul.to_tree(df)

    return final_tree
    # return cnt_tree1, cnt_tree2
