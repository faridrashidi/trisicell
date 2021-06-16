from typing import Optional

import metools as met
import networkx as nx
import numpy as np


def fitch(tree: nx.DiGraph, seed: Optional[int] = None):
    """[summary].

    Parameters
    ----------
    tree : nx.DiGraph
        The tree in which leaves has `'O'` attribute.
    seed : Optional[int], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """
    tree_bu = _fitch_hartigan_bottom_up(tree)
    tree_td = _fitch_hartigan_top_down(tree_bu, seed=seed)
    tree_fi = _parsimony_cost(tree_td)
    return tree_fi, tree_td, tree_bu


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


def _fill_child(a, b):
    if b == 2:
        return a
    return b


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


_get_consensus = np.vectorize(_get_consensus)
_fill_child = np.vectorize(_fill_child)


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
            gc = _get_consensus(
                tree.nodes[children[0]]["V"], tree.nodes[children[1]]["V"]
            )
            tree.nodes[i]["V"] = gc
            # if len(children) > 2:
            #     print(children)
            #     print(tree.nodes[children[0]]["V"], tree.nodes[children[1]]["V"])
            #     print(gc)
        d -= 1

    return tree


def _set_assignment(v, seed):
    if seed is not None:
        np.random.seed(seed)
    idx = np.where(v == 2)[0]
    for i in idx:
        v[i] = np.random.choice([0, 1])
    return v


def _fitch_hartigan_top_down(tree2, seed):
    tree = tree2.copy()
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
    # leaves = [n for n in tree.nodes() if tree.out_degree(n) == 0]

    tree = _set_depth(tree, root)
    md = _get_max_depth(tree, root)
    tree.nodes[root]["V"] = _set_assignment(tree.nodes[root]["V"], seed=seed)
    d = 1

    while d <= md:
        internal_nodes = list(_cut_tree(tree, d))
        for i in internal_nodes:
            parent = list(tree.predecessors(i))[0]
            tree.nodes[i]["V"] = _fill_child(
                tree.nodes[parent]["V"], tree.nodes[i]["V"]
            )
        d += 1
    return tree


def _parsimony_cost(tree):
    cost = 0
    for e in tree.edges():
        source = e[0]
        dest = e[1]
        s = np.sum(np.abs(tree.nodes[source]["V"] - tree.nodes[dest]["V"]))
        if s > 0:
            cost += s
            tree.edges[e]["F"] = np.where(
                tree.nodes[source]["V"] != tree.nodes[dest]["V"]
            )[0]
    met.logg.info(f"cost: {cost}")
    return tree
