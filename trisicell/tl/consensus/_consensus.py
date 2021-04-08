import itertools

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc


def consensus_combine(df):
    """Combine the replicates/cells in genotype matrix.

    This function combine the genotype of replicates/cells that are
    in the form of `{CellA}_{id}` to `{CellA}`.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        The input genotype matrix in conflict-free format.

    Returns
    -------
    :class:`pandas.DataFrame`
        The combine genotype matrix.
    """

    df2 = df.groupby(df.index.str.split("_").str[0]).transform("prod")
    df2 = df2.groupby(df2.index.str.split("_").str[0]).first()
    return df2


def _consensus_of_cfmatrices(df1, df2):
    # TODO: unused
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


def _get_cnt_tree(tree):
    cnt_tree = tree.copy()
    for u, v, l in cnt_tree.edges.data("label"):
        mutations = l.split(cnt_tree.graph["splitter_mut"])
        cnt_tree.add_edge(u, v, label=len(mutations), mutations=mutations)
        if cnt_tree.in_degree(v) != 0:
            if "––" not in tree.nodes[v]["label"]:
                cnt_tree.nodes[v]["label"] = cnt_tree.nodes[v]["label"].split(
                    cnt_tree.graph["splitter_cell"]
                )
    return cnt_tree


def _get_node_labels(tree):
    nl = {}
    for x in tree.nodes:
        if tree.in_degree(x) != 0:
            if "––" not in tree.nodes[x]["label"]:
                for y in tree.nodes[x]["label"]:
                    nl[y] = x
        else:
            nl["root"] = x
    return nl


def _has_path(tree, nodes, l1, l2):
    try:
        n1 = nodes[l1]
        n2 = nodes[l2]
        dij = nx.bidirectional_dijkstra(tree, n1, n2, weight="label")
        return True, dij[0]
    except:
        return False, None
    return False, None


def _merge_with_parents(tree, nodes, n1, n2):
    lca = nx.lowest_common_ancestor(tree, n1, n2)
    costs, paths = nx.bidirectional_dijkstra(tree, lca, n1, weight="label")
    paths = paths[::-1]
    labels = []
    for p in paths[1:]:
        #### merge-with-parent, n1 -> n2 (n2 is going to be replaced by n1)
        if "––" not in tree.nodes[p]["label"]:
            if "––" not in tree.nodes[n1]["label"]:
                tree.nodes[n1]["label"] += tree.nodes[p]["label"]
            else:
                tree.nodes[n1]["label"] = tree.nodes[p]["label"]
            for x in tree.nodes[p]["label"]:
                nodes[x] = n1
        cost = tree.edges[p, n1]["label"]

        #### add the mutations of this edge to the grandparent edge
        grandparent = list(tree.predecessors(p))[0]
        tree.edges[(grandparent, p)]["mutations"] += tree.edges[(p, n1)]["mutations"]
        tree.edges[(grandparent, p)]["label"] += tree.edges[(p, n1)]["label"]

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
                nx.algorithms.traversal.depth_first_search.dfs_tree(cnt_tree, x).nodes,
            )
            for y in sub.nodes:
                if "––" not in cnt_tree.nodes[y]["label"]:
                    lbs += cnt_tree.nodes[y]["label"]
            labels[x] = lbs
    return labels


def _add_private_muts(cnt_tree, sc_data, tree_nodes):
    latest = max(cnt_tree.nodes)
    nodes = list(cnt_tree.nodes)
    has_expanded_all = []
    for node in nodes:
        if "––" not in cnt_tree.nodes[node]["label"] and cnt_tree.in_degree(node) != 0:
            has_expanded = []
            for cell in cnt_tree.nodes[node]["label"]:
                private = (sc_data.loc[cell] == 1) & (
                    sc_data.loc[np.setdiff1d(sc_data.index, cell)].sum() == 0
                )
                private = list(private[private == True].index)
                if len(private) > 0:
                    has_expanded_all.append(cell)
                    has_expanded.append(cell)
                    latest = latest + 1
                    cnt_tree.add_node(latest, label=[cell])
                    cnt_tree.add_edge(
                        node, latest, label=len(private), mutations=private
                    )
                    tree_nodes[cell] = latest
            for cell in has_expanded:
                cnt_tree.nodes[node]["label"].remove(cell)
            if len(cnt_tree.nodes[node]["label"]) == 0:
                cnt_tree.nodes[node]["label"] = "––"
    return cnt_tree, tree_nodes, has_expanded_all


def consensus_tree(sc1, sc2):
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

    #### remove private mutations and build the trees
    sc1_c = sc1.loc[:, sc1.sum() > 1].copy()
    sc2_c = sc2.loc[:, sc2.sum() > 1].copy()

    tree1 = tsc.ul.to_tree(sc1_c)
    tree2 = tsc.ul.to_tree(sc2_c)

    cnt_tree1 = _get_cnt_tree(tree1)
    cnt_tree2 = _get_cnt_tree(tree2)

    nodes1 = _get_node_labels(cnt_tree1)
    nodes2 = _get_node_labels(cnt_tree2)

    common_cells = list(np.intersect1d(nodes1.keys(), nodes2.keys())[0])

    total_cost = 0

    #### step 1,2
    for x, y in itertools.permutations(common_cells, 2):
        hp1, cost1 = _has_path(cnt_tree1, nodes1, x, y)
        if hp1:
            hp2, cost2 = _has_path(cnt_tree2, nodes2, x, y)
            if not hp2:
                cnt_tree2, nodes2, cost2 = _merge_with_parents(
                    cnt_tree2, nodes2, nodes2[x], nodes2[y]
                )
                tsc.logg.info(f"merge in tree2: `{x}` to `{y}` and cost={cost2}")
                total_cost += cost2

    for x, y in itertools.permutations(common_cells, 2):
        hp2, cost2 = _has_path(cnt_tree2, nodes2, x, y)
        if hp2:
            hp1, cost1 = _has_path(cnt_tree1, nodes1, x, y)
            if not hp1:
                cnt_tree1, nodes1, cost1 = _merge_with_parents(
                    cnt_tree1, nodes1, nodes1[x], nodes1[y]
                )
                tsc.logg.info(f"merge in tree1: `{x}` to `{y}` and cost={cost1}")
                total_cost += cost1

    for x, y in itertools.permutations(common_cells, 2):
        hp1, cost1 = _has_path(cnt_tree1, nodes1, x, y)
        if hp1:
            hp2, cost2 = _has_path(cnt_tree2, nodes2, x, y)
            if not hp2:
                cnt_tree2, nodes2, cost2 = _merge_with_parents(
                    cnt_tree2, nodes2, nodes2[x], nodes2[y]
                )
                tsc.logg.info(f"merge in tree2: `{x}` to `{y}` and cost={cost2}")
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
        cnt_tree1, nodes1, cost1 = _merge_with_parents(cnt_tree1, nodes1, y, parent)
        tsc.logg.info(
            f"merge in tree1 - internal: `node{y}` to `node{parent}` and cost={cost1}"
        )
        total_cost += cost1

    for y, _ in labels2.items():
        parent = list(cnt_tree2.predecessors(y))[0]
        cnt_tree2, nodes2, cost2 = _merge_with_parents(cnt_tree2, nodes2, y, parent)
        tsc.logg.info(
            f"merge in tree2 - internal: `node{y}` to `node{parent}` and cost={cost2}"
        )
        total_cost += cost2

    #### add private mutations back
    cnt_tree1, nodes1, has_expanded1 = _add_private_muts(cnt_tree1, sc1, nodes1)
    cnt_tree2, nodes2, has_expanded2 = _add_private_muts(cnt_tree2, sc2, nodes2)

    for y in np.setdiff1d(has_expanded1, has_expanded2):
        parent = list(cnt_tree1.predecessors(nodes1[y]))[0]
        cnt_tree1, nodes1, cost1 = _merge_with_parents(
            cnt_tree1, nodes1, nodes1[y], parent
        )
        tsc.logg.info(
            f"merge in tree1 - private: `{y}` to `node{parent}` and cost={cost1}"
        )
        total_cost += cost1

    for y in np.setdiff1d(has_expanded2, has_expanded1):
        parent = list(cnt_tree2.predecessors(nodes2[y]))[0]
        cnt_tree2, nodes2, cost2 = _merge_with_parents(
            cnt_tree2, nodes2, nodes2[y], parent
        )
        tsc.logg.info(
            f"merge in tree2 - private: `{y}` to `node{parent}` and cost={cost2}"
        )
        total_cost += cost1

    tsc.logg.info("    ---> total cost:", total_cost)

    if nx.is_isomorphic(cnt_tree1, cnt_tree2) == False:
        raise RuntimeError("Error: Two trees are not isomorphic!")

    final_tree = cnt_tree1.copy()
    root1 = [x for x in cnt_tree1.nodes if cnt_tree1.in_degree(x) == 0][0]
    root2 = [x for x in cnt_tree2.nodes if cnt_tree2.in_degree(x) == 0][0]
    leaves2 = [x for x in cnt_tree2.nodes if cnt_tree2.out_degree(x) == 0]
    for leaf2 in leaves2:
        list2 = nx.dijkstra_path(cnt_tree2, root2, leaf2)
        label1 = cnt_tree2.nodes[leaf2]["label"]

        mutations2 = []
        for first, second in zip(list2, list2[1:]):
            mutations2.append(cnt_tree2.edges[(first, second)]["mutations"])

        list1 = nx.dijkstra_path(cnt_tree1, root1, nodes1[label1[0]])
        i = 0
        for first, second in zip(list1, list1[1:]):
            mutations1 = cnt_tree1.edges[(first, second)]["mutations"]
            mutations = np.intersect1d(mutations1, mutations2[i])
            final_tree.add_edge(
                first, second, label=final_tree.graph["splitter_mut"].join(mutations)
            )
            i += 1
    for v in final_tree.nodes:
        if final_tree.in_degree(v) != 0:
            if "––" not in final_tree.nodes[v]["label"]:
                final_tree.nodes[v]["label"] = final_tree.graph["splitter_cell"].join(
                    final_tree.nodes[v]["label"]
                )

    data = tsc.ul.to_cfmatrix(final_tree)
    final_tree.graph["data"] = data
    del final_tree.graph["become_germline"]

    return final_tree
