from collections import Counter, defaultdict
from functools import partial
from itertools import combinations
from multiprocessing import Pool, set_start_method

import networkx as nx
import numpy as np

__author__ = "Simone Ciccolella"
__date__ = "9/24/21"


class Tree:
    def __init__(self, T, label_to_nodes, node_to_labels, label_set):
        self.T = T
        self.label_to_nodes = label_to_nodes
        self.node_to_labels = node_to_labels
        self.label_set = label_set

        self.LCA = LCA(self.T, self.label_to_nodes, self.node_to_labels)


class LCA:
    def __init__(self, T, T_label_to_node, T_node_to_labels):
        self.LCA_dict = {}
        self.LCA_label_dict = T_label_to_node
        self.LCA_node2lbl = T_node_to_labels
        for lca in nx.tree_all_pairs_lowest_common_ancestor(T):
            self.LCA_dict[(lca[0][0], lca[0][1])] = lca[1]

    def lca_nodes(self, node1, node2):
        try:
            return self.LCA_dict[node1, node2]
        except Exception:
            return self.LCA_dict[node2, node1]

    def lca_labels(self, label1, label2):
        nodes1 = self.LCA_label_dict[label1]
        nodes2 = self.LCA_label_dict[label2]

        lca_mset = Counter()
        for n1 in nodes1:
            for n2 in nodes2:
                lca_mset.update(self.lca_nodes(n1, n2))

        return lca_mset

    def label_to_node(self, label):
        return self.LCA_label_dict[label]

    def node_to_labels(self, node):
        return self.LCA_node2lbl[node]

    def __str__(self):
        """Return."""
        return str(self.LCA_dict)


class ExtValue:
    def __init__(self):
        self.counter = 1
        self.nodes = defaultdict(set)

    def get(self, node):
        if node not in self.nodes:
            self.nodes[node].add(f"EXT{self.counter}")
            self.counter += 1
        return self.nodes[node]


def sigmoid(x, mult=10.0):
    if x == 0:
        return 0
    if x == 1:
        return 1
    return 1 / (1 + np.exp(-mult * (x - 0.5)))


def intersect_mset_card(list_lca1, list_lca2):
    mset1 = defaultdict(int)
    mset2 = defaultdict(int)

    for lca in list_lca1:
        mset1[",".join(str(x) for x in lca)] += 1

    for lca in list_lca2:
        mset2[",".join(str(x) for x in lca)] += 1

    card = 0

    for k in mset1:
        if k in mset2:
            card += min(mset1[k], mset2[k])

    return card


def is_equal_struct(triple, LCA1, LCA2):
    t = sorted(triple)
    t_set = set(triple)

    triples_nodes_T1 = []

    for node1 in LCA1.label_to_node(t[0]):
        for node2 in LCA1.label_to_node(t[1]):
            for node3 in LCA1.label_to_node(t[2]):
                triples_nodes_T1.append([node1, node2, node3])

    triples_nodes_T2 = []

    for node1 in LCA2.label_to_node(t[0]):
        for node2 in LCA2.label_to_node(t[1]):
            for node3 in LCA2.label_to_node(t[2]):
                triples_nodes_T2.append([node1, node2, node3])

    cmp_vecs_t1 = [[None, None, None] for x in range(len(triples_nodes_T1))]
    ext_t1 = [ExtValue() for x in range(len(triples_nodes_T1))]

    cmp_vecs_t2 = [[None, None, None] for x in range(len(triples_nodes_T2))]
    ext_t2 = [ExtValue() for x in range(len(triples_nodes_T2))]

    for ix_c, couple in enumerate(combinations(t, 2)):
        ix_lb1 = t.index(couple[0])
        ix_lb2 = t.index(couple[1])
        ix_other = list({0, 1, 2} - {ix_lb1, ix_lb2})[0]

        for ix_t1, triple_T1 in enumerate(triples_nodes_T1):
            lca_nd = LCA1.lca_nodes(triple_T1[ix_lb1], triple_T1[ix_lb2])
            lca_lb = LCA1.node_to_labels(lca_nd)
            if len(lca_lb & t_set) > 0:
                cmp_vecs_t1[ix_t1][ix_c] = lca_lb & t_set
            else:
                lca_nd_other = LCA1.lca_nodes(lca_nd, triple_T1[ix_other])
                lca_lb_other = LCA1.node_to_labels(lca_nd_other)
                if len(lca_lb_other & t_set) > 0:
                    cmp_vecs_t1[ix_t1][ix_c] = lca_lb_other & t_set
                else:
                    cmp_vecs_t1[ix_t1][ix_c] = ext_t1[ix_t1].get(lca_nd)

        for ix_t2, triple_T2 in enumerate(triples_nodes_T2):
            lca_nd = LCA2.lca_nodes(triple_T2[ix_lb1], triple_T2[ix_lb2])
            lca_lb = LCA2.node_to_labels(lca_nd)
            if len(lca_lb & t_set) > 0:
                cmp_vecs_t2[ix_t2][ix_c] = lca_lb & t_set
            else:
                lca_nd_other = LCA2.lca_nodes(lca_nd, triple_T2[ix_other])
                lca_lb_other = LCA2.node_to_labels(lca_nd_other)
                if len(lca_lb_other & t_set) > 0:
                    cmp_vecs_t2[ix_t2][ix_c] = lca_lb_other & t_set
                else:
                    cmp_vecs_t2[ix_t2][ix_c] = ext_t2[ix_t2].get(lca_nd)

    missing = True if len(cmp_vecs_t1) == 0 or len(cmp_vecs_t2) == 0 else False

    return (
        missing,
        intersect_mset_card(cmp_vecs_t1, cmp_vecs_t2),
        max(len(cmp_vecs_t1), len(cmp_vecs_t2)),
    )


def get_nset_sig(x_i, x_u):
    return x_u + sigmoid(x_i) * (x_i - x_u)


def similarity(tree1, tree2, mode="sigmoid", sigmoid_mult=10.0, cores=1):
    """
    Compute the similarity score of the two trees.

    Parameters:
    tree1 (Tree): MP3-treesim tree representation
    tree2 (Tree): MP3-treesim tree representation

    Keyword arguments:
    mode (str): 'sigmoid', 'intersection', 'union'
                 or 'geometric',
                 sets the similarity calculation.
                 By default is set to 'sigmoid'.
    sigmoid_mult (float): Multiplicator for the
                 sigmoid calculation.
    cores (int); Number of cores used for the computation.

    Returns:
    float: Similarity score
    """

    if mode not in ["sigmoid", "intersection", "union", "geometric"]:
        raise AttributeError("Incorrect value of mode passed.")

    numerator = 0
    denominator_i = 0
    denominator_u = 0

    if len(set(tree1.label_set) & set(tree2.label_set)) == 0:
        return 0.0

    if mode == "intersection":
        labels = set(tree1.label_set) & set(tree2.label_set)
    else:
        labels = set(tree1.label_set) | set(tree2.label_set)

    combs = combinations(labels, 3)

    if cores <= 0:
        cores = None

    if cores == 1:
        for triple in combs:
            missing, num, dem = is_equal_struct(triple, tree1.LCA, tree2.LCA)
            numerator += num
            if missing:
                denominator_u += dem
            else:
                denominator_i += dem
    else:
        try:
            set_start_method("spawn")
        except RuntimeError:
            pass

        with Pool(processes=cores) as pool:
            func = partial(is_equal_struct, LCA1=tree1.LCA, LCA2=tree2.LCA)
            results = pool.map(func, combs)

            for res in results:
                missing, num, dem = res
                numerator += num
                if missing:
                    denominator_u += dem
                else:
                    denominator_i += dem

    similarity_score = 0
    if mode == "intersection":
        similarity_score = float(numerator) / denominator_i
    elif mode == "union":
        similarity_score = float(numerator) / (denominator_i + denominator_u)
    else:
        similarity_score_i = float(numerator) / denominator_i
        similarity_score_u = float(numerator) / (denominator_i + denominator_u)

        if mode == "sigmoid":
            similarity_score = similarity_score_u + sigmoid(
                similarity_score_i, mult=sigmoid_mult
            ) * min(similarity_score_i - similarity_score_u, similarity_score_u)

        elif mode == "geometric":
            similarity_score = np.sqrt(similarity_score_i * similarity_score_u)

    return similarity_score


def build_tree(T, labeled_only=False, exclude=None):
    """
    Build the MP3-treesim tree representation from a networkx representation.

    The tree must have a attribute `label` for each node. Labels in a node
    must be separated by a comma.

    Parameters:
    T (nx.nx_agraph): Tree in networkx representation
    labeled_only (bool): If true nodes without attribute `label`
    will be ignored, meaning that T is a partially labeled tree.
    exclude (list(str)): List of labels to exclude from computation

    Returns:
    Tree: MP3-treesim tree representation
    """

    if not nx.is_tree(T):
        raise ValueError("Not a valid tree.")

    label_to_nodes = defaultdict(set)
    label_set = set()
    node_to_labels = defaultdict(set)

    for node in T.nodes(data=True):
        id_node = node[0]

        if "label" not in node[1] and not labeled_only:
            node[1]["label"] = str(node[0])
        if "label" not in node[1] and labeled_only:
            node[1]["label"] = ""
            continue

        labels = node[1]["label"].split(",")
        for label in labels:
            if exclude:
                if label in exclude:
                    continue

            label_set.add(label)
            label_to_nodes[label].add(id_node)
            node_to_labels[id_node].add(label)

    label_set = list(label_set)
    label_set = sorted(label_set)

    return Tree(T, label_to_nodes, node_to_labels, label_set)
