from operator import index
from sys import path

import networkx as nx
import numpy as np
import pandas as pd

import trisicell as tsc


def to_tree(df):
    """Convert a conflict-free matrix to a tree object.

    This function converts a conflict-free matrix to a tree object in which
    nodes are labled with cells and edges are lables with mutations. The root is
    labled by 'root'. Mutations are seperated by `.graph['splitter_mut']` and cells
    are seperated by `.graph['splitter_cell']`. Those mutations that are not present
    in any cell are stored in `.graph['become_germline']`. Mutations happed once
    during the evolution so there is no repetitive mutation.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        A genotype dataframe in which rows are cells and columns are mutations.
        Note that this dataframe must be conflict-free.

    Returns
    -------
    :class:`networkx.DiGraph`
        A perfect phylogenetic tree.

    Raises
    ------
    Exception
        If the input dataframe is not conflict-free.
    """

    if not tsc.ul.is_conflict_free_gusfield(df):
        raise Exception("The input is not conflict-free!")

    def _contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    tree = nx.DiGraph()
    tree.graph["data"] = df
    tree.graph["splitter_mut"] = "\n"
    tree.graph["splitter_cell"] = "\n"
    tree.graph["become_germline"] = df.columns[(df == 0).all(axis=0)]

    matrix = df.values
    names_mut = list(df.columns)

    i = 0
    while i < matrix.shape[1]:
        j = i + 1
        while j < matrix.shape[1]:
            if np.array_equal(matrix[:, i], matrix[:, j]):
                matrix = np.delete(matrix, j, 1)
                x = names_mut.pop(j)
                names_mut[i] += tree.graph["splitter_mut"] + x
                j -= 1
            j += 1
        i += 1

    rows = matrix.shape[0]
    cols = matrix.shape[1]
    dimensions = np.sum(matrix, axis=0)
    indices = np.argsort(dimensions)
    dimensions = np.sort(dimensions)
    names_mut = [names_mut[indices[i]] for i in range(cols)]

    tree.add_node(cols)
    tree.add_node(cols - 1)
    tree.add_edge(cols, cols - 1, label=names_mut[cols - 1])
    node_mud = {}
    node_mud[names_mut[cols - 1]] = cols - 1

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break
        attached = False
        for j in range(i + 1, cols):
            if _contains(matrix[:, indices[j]], matrix[:, indices[i]]):
                tree.add_node(i)
                tree.add_edge(node_mud[names_mut[j]], i, label=names_mut[i])
                node_mud[names_mut[i]] = i
                attached = True
                break
        if not attached:
            tree.add_node(i)
            tree.add_edge(cols, i, label=names_mut[i])
            node_mud[names_mut[i]] = i
        i -= 1

    tumor_cells = []
    clusters = {cols: "root"}
    for node in tree:
        if node == cols:
            tree.nodes[node]["label"] = "root"
            continue
        untilnow_mut = []
        sp = nx.shortest_path(tree, cols, node)
        for i in range(len(sp) - 1):
            untilnow_mut += tree.get_edge_data(sp[i], sp[i + 1])["label"].split(
                tree.graph["splitter_mut"]
            )
        untilnow_cell = df.loc[
            (df[untilnow_mut] == 1).all(axis=1)
            & (df[[x for x in df.columns if x not in untilnow_mut]] == 0).all(axis=1)
        ].index
        if len(untilnow_cell) > 0:
            clusters[node] = f"{tree.graph['splitter_cell'].join(untilnow_cell)}"
            tumor_cells += [y for y in tree.graph["splitter_cell"].join(untilnow_cell)]
        else:
            clusters[node] = "––"

        tree.nodes[node]["label"] = clusters[node]

    tree.graph["normal_cells"] = df[df.sum(axis=1) == 0].index
    tree.graph["root_id"] = cols

    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1
    return tree


def to_cfmatrix(tree):
    """Convert phylogenetic tree to conflict-free matrix.

    This function converts a phylogenetic tree in which mutations
    are at edges and cells are at nodes to a conflict-free matrix
    where rows are cells, columns are mutations and each entry is
    either zero or one representing the absence or presence of the
    mutation in the cell.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The phylogenetic tree.

    Returns
    -------
    :class:`pandas.DataFrame`
        The conflict-free matrix.
    """

    mutations = []
    cells = []
    for _, v, l in tree.edges(data=True):
        mutations += l["label"].split(tree.graph["splitter_mut"])
        if "––" not in tree.nodes[v]["label"]:
            cells += tree.nodes[v]["label"].split(tree.graph["splitter_cell"])
    df = pd.DataFrame(0, index=cells, columns=mutations)
    root = tsc.ul.root_id(tree)
    leaves = [x for x in tree.nodes if tree.out_degree(x) == 0]
    for leaf in leaves:
        nodes = nx.dijkstra_path(tree, root, leaf)
        mut = []
        for first, second in zip(nodes, nodes[1:]):
            mut += tree.edges[(first, second)]["label"].split(
                tree.graph["splitter_mut"]
            )
            if "––" not in tree.nodes[second]["label"]:
                cell = tree.nodes[second]["label"].split(tree.graph["splitter_cell"])
                df.loc[cell, mut] = 1
    return df


def to_mtree(tree):
    """Convert the phylogenetic tree to mutation tree.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The phylogenetic tree in which cells are in nodes and
        mutations are at edges.

    Returns
    -------
    :class:`networkx.DiGraph`
        The mutation tree in which mutations are in nodes.
    """

    tree2 = nx.DiGraph()
    for u, v, l in tree.edges.data("label"):
        if tree.in_degree(u) == 0:
            tree2.add_node(u, label="root")
        muts = l.split(tree.graph["splitter_mut"])
        tree2.add_node(v, label=muts)
        tree2.add_edge(u, v)
    return tree2


def _to_newick(tree):
    def subtree(at):
        return nx.subgraph(
            tree,
            nx.algorithms.traversal.depth_first_search.dfs_tree(tree, at).nodes - [at],
        )

    def children(at):
        return [n for n in tree.neighbors(at)]

    root = tsc.ul.root_id(tree)

    def newick_recursive(node_id):
        node_ids = children(node_id)
        if len(node_ids) == 0:
            cells = tree.nodes[node_id]["label"].split(tree.graph["splitter_cell"])
            return "(" + ":1,".join(cells) + f":1)Node{node_id+1}:1"
        elif len(node_ids) > 0:
            cells = tree.nodes[node_id]["label"].split(tree.graph["splitter_cell"])
            if not ("––" in cells or "root" in cells):
                return (
                    "("
                    + ":1,".join(cells)
                    + ":1,"
                    + ",".join([newick_recursive(node_id) for node_id in node_ids])
                    + f":1)Node{node_id+1}"
                )
            else:
                return (
                    "("
                    + ":1,".join([newick_recursive(node_id) for node_id in node_ids])
                    + f":1)Node{node_id+1}:1"
                )

    newick = newick_recursive(root) + ";"
    return newick


def _split_labels(mt, mt_guide):
    root_guide = [node for node in mt_guide.nodes if mt_guide.in_degree(node) == 0][0]
    guide = [
        y
        for x in nx.algorithms.bfs_tree(mt_guide, root_guide).nodes
        for y in mt_guide.nodes[x]["label"]
    ]

    root = [node for node in mt.nodes if mt.in_degree(node) == 0][0]
    nodes = nx.algorithms.bfs_tree(mt, root).nodes

    latest_node = root
    removing = []
    for node in nodes:
        muts = mt.nodes[node]["label"]
        if len(muts) > 1 and node != root:
            parent = list(mt.predecessors(node))[-1]
            if mt.out_degree(node) != 0:
                children = list(mt.successors(node))
            else:
                children = -1
            muts = sorted(muts, key=lambda i: guide.index(i))

            removing.append(node)
            for i, x in enumerate(muts):
                latest_node += 1
                mt.add_node(latest_node, label=x)
                if i == 0:
                    mt.add_edge(parent, latest_node)
                elif i == len(muts) - 1:
                    if children != -1:
                        mt.add_edge(latest_node - 1, latest_node)
                        for child in children:
                            mt.add_edge(latest_node, child)
                    else:
                        mt.add_edge(latest_node - 1, latest_node)
                else:
                    mt.add_edge(latest_node - 1, latest_node)
        else:
            if node != root:
                mt.nodes[node]["label"] = mt.nodes[node]["label"][0]

    for node in removing:
        mt.remove_node(node)
    return mt


def _to_apted(sl_tree):
    def children(at):
        return [n for n in sl_tree.neighbors(at)]

    def apted_recursive(node):
        nodes = children(node)
        if len(nodes) == 0:
            l = sl_tree.nodes[node]["label"]
            return "{" + l + "}"
        else:
            l = sl_tree.nodes[node]["label"]
            x = ""
            for node in nodes:
                x += apted_recursive(node)
            return "{" + l + x + "}"

    root = [node for node in sl_tree.nodes if sl_tree.in_degree(node) == 0][0]
    return apted_recursive(root)


def root_id(tree):
    return [x for x in tree.nodes if tree.in_degree(x) == 0][0]


# def partition_cells(tree, node):
#     nd = tree.graph["splitter_cell"].join(node)
#     cells = []
#     for x in list(nx.algorithms.traversal.depth_first_search.dfs_tree(tree, nd).nodes):
#         for y in x.split(", "):
#             if not "–" in y:
#                 cells.append(y)
#     cells = np.array(cells)
#     return cells, np.setdiff1d(tree.graph["data"].index, cells)


def cells_rooted_at(tree, node_id):
    muts = tree.graph["mutation_list"][tree.graph["mutation_list"].Node == node_id]
    if muts.index.shape[0] == 0:
        return np.array([])
    cells = (tree.graph["data"][muts.index] == 1).all(axis=1)
    cells = np.array(tree.graph["data"].loc[cells].index)
    return cells  # , np.setdiff1d(tree.graph["data"].index, cells)


def muts_rooted_at(tree, node_id):
    muts = tree.graph["mutation_list"][tree.graph["mutation_list"].Node == node_id]
    if muts.index.shape[0] == 0:
        return np.array([])

    muts = muts.index.tolist()
    nd = int(node_id.replace("[", "").replace("]", ""))
    paths = nx.algorithms.traversal.depth_first_search.dfs_tree(tree, nd).nodes
    sub_tree = nx.subgraph(tree, paths)
    for u, v, l in sub_tree.edges.data("label"):
        muts += l.split(tree.graph["splitter_mut"])
    return np.array(muts)
