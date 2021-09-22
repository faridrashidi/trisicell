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
    last_node = max(nodes)
    for n in nodes:
        if cp_tree.in_degree(n) == 1 and cp_tree.out_degree(n) > 0:
            if "––" not in cp_tree.nodes[n]["label"]:
                for x in cp_tree.nodes[n]["label"]:
                    artificial_added.append(x)
                    last_node += 1
                    cp_tree.add_node(last_node, label=[x])
                    cp_tree.add_edge(n, last_node)
                cp_tree.nodes[n]["label"] = "––"
        elif cp_tree.in_degree(n) == 1 and cp_tree.out_degree(n) == 0:
            if len(cp_tree.nodes[n]["label"]) > 1:
                for x in cp_tree.nodes[n]["label"]:
                    artificial_added.append(x)
                    last_node += 1
                    cp_tree.add_node(last_node, label=[x])
                    cp_tree.add_edge(n, last_node)
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

    # step 1
    cnt_tree1, artificial_added1 = _add_artificial_edge(cnt_tree1)
    cnt_tree2, artificial_added2 = _add_artificial_edge(cnt_tree2)
    artificial_added = list(np.union1d(artificial_added1, artificial_added2))
    tsc.logg.debug(artificial_added)

    # step 2
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

    # step 3
    final_tree = _get_cnt_tree(final_tree)
    final_tree, _ = _add_artificial_edge(final_tree)

    # step 4
    nodes = list(final_tree.nodes)
    for n in nodes:
        if "––" not in final_tree.nodes[n]["label"] and final_tree.in_degree(n) == 1:
            if final_tree.nodes[n]["label"][0] in artificial_added:
                m = list(final_tree.predecessors(n))[0]
                if "––" in final_tree.nodes[m]["label"]:
                    final_tree.nodes[m]["label"] = final_tree.nodes[n]["label"]
                else:
                    final_tree.nodes[m]["label"] += final_tree.nodes[n]["label"]
                final_tree = nx.contracted_nodes(final_tree, m, n, self_loops=False)

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
