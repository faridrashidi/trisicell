import networkx as nx
import numpy as np

import trisicell as tsc


def fitch(tree: nx.DiGraph):
    """[summary].

    Parameters
    ----------
    tree : nx.DiGraph
        The tree in which leaves has `'profile'` attribute.

    Returns
    -------
    [type]
        [description]
    """
    tree_bu = _fitch_hartigan_bottom_up(tree)
    tree_td = _fitch_hartigan_top_down(tree_bu)
    tree_fi = _parsimony_cost(tree_td)
    return tree_fi


def _get_consensus(a, b):
    if a == 1 and b == 1:
        return 1
    if a == 0 and b == 0:
        return 0
    if a == 2 and (b == 0 or b == 1):
        return b
    if b == 2 and (a == 0 or a == 1):
        return a
    return 2


_get_consensus = np.vectorize(_get_consensus)


def _fill_child(a, b):
    if b == 2:
        return a
    return b


_fill_child = np.vectorize(_fill_child)


def _set_depth(G, root):
    depth = nx.shortest_path_length(G, root)
    for d in depth.keys():
        G.nodes[d]["depth"] = depth[d]
    return G


def _cut_tree(G, depth):
    nodes = []
    for n in G.nodes:
        if G.nodes[n]["depth"] == depth:
            nodes.append(n)
    return nodes


def _get_max_depth(G, root):
    md = 0
    for n in nx.descendants(G, root):
        if G.nodes[n]["depth"] > md:
            md = G.nodes[n]["depth"]
    return md


def _fitch_hartigan_bottom_up(tree2):
    tree = tree2.copy()
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
    # leaves = [n for n in tree.nodes() if tree.out_degree(n) == 0]

    tree = _set_depth(tree, root)
    md = _get_max_depth(tree, root)
    d = md

    while d >= 0:
        internal_nodes = _cut_tree(tree, d)
        for i in internal_nodes:
            children = list(tree.successors(i))
            if len(children) == 0:
                continue
            if len(children) == 1:
                tree.nodes[i]["profile"] = tree.nodes[children[0]]["profile"]
            else:
                gc = tree.nodes[children[0]]["profile"]
                for child in children[1:]:
                    gc = _get_consensus(
                        gc,
                        tree.nodes[child]["profile"],
                    )
                tree.nodes[i]["profile"] = gc
        d -= 1

    return tree


def _fitch_hartigan_top_down(tree2):
    tree = tree2.copy()
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]

    tree = _set_depth(tree, root)
    md = _get_max_depth(tree, root)
    tree.nodes[root]["profile"] = [0] * tree.nodes[root]["profile"].shape[0]
    d = 1

    while d <= md:
        internal_nodes = list(_cut_tree(tree, d))
        for i in internal_nodes:
            parent = list(tree.predecessors(i))[0]
            tree.nodes[i]["profile"] = _fill_child(
                tree.nodes[parent]["profile"], tree.nodes[i]["profile"]
            )
        d += 1
    return tree


def _parsimony_cost(tree):
    cost = 0
    for e in tree.edges():
        source = e[0]
        dest = e[1]
        s = np.sum(np.abs(tree.nodes[source]["profile"] - tree.nodes[dest]["profile"]))
        if s > 0:
            cost += s
            tree.edges[e]["event"] = np.where(
                tree.nodes[source]["profile"] != tree.nodes[dest]["profile"]
            )[0]
        else:
            tree.edges[e]["event"] = []
    tsc.logg.info(f"cost: {cost}")
    return tree
