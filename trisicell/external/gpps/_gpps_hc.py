import random
from functools import lru_cache

import numpy as np

import trisicell as tsc
from trisicell.external.gpps._utils_hc import (
    copy_tree,
    import_ilp_out,
    prune_and_reattach,
)

__author__ = "Simone Ciccolella"
__date__ = "11/30/21"


@lru_cache(maxsize=None)
def cell_row_likelihood(input_row, node_genotype, alpha, beta):
    likelihood = 0
    for j in range(len(input_row)):
        if input_row[j] == "0":
            if node_genotype[j] == "0":
                likelihood += np.log(1 - beta)
            elif node_genotype[j] == "1":
                likelihood += np.log(alpha)
            else:
                likelihood += -np.inf
        elif input_row[j] == "1":
            if node_genotype[j] == "0":
                likelihood += np.log(beta)
            elif node_genotype[j] == "1":
                likelihood += np.log(1 - alpha)
            else:
                likelihood += -np.inf

        elif input_row[j] == "2":
            likelihood += 0
    return likelihood


def greedy_tree_likelihood(tree, nid_dict, input_scs, alpha, beta):
    tree
    likelihood = 0
    attachment = []
    for row in input_scs:
        str_row = "".join(str(int(x)) for x in row)
        best_lh = -np.inf
        best_attachment = -1
        for node in nid_dict:
            str_gt = "".join(str(int(x)) for x in nid_dict[node].genotype_profile)

            lh = cell_row_likelihood(str_row, str_gt, alpha, beta)
            if lh > best_lh:
                best_lh = lh
                best_attachment = node
        likelihood += best_lh
        attachment.append(best_attachment)

    return likelihood, attachment


def get_expect_matrix(tree, nid_dict, input_scs, alpha, beta):
    _, attachment = greedy_tree_likelihood(tree, nid_dict, input_scs, alpha, beta)
    e_matrix = []

    for node in nid_dict:
        t_node = nid_dict[node]
        if (
            t_node.loss
            and t_node.id_node not in attachment
            and len(t_node.children) == 0
        ):
            tsc.logg.debug(t_node.name)

    for c in range(len(attachment)):
        e_matrix.append(nid_dict[attachment[c]].genotype_profile)
    return e_matrix


def generate_neighborhood(start_tree, start_nid_dict, neighborhood_size):
    start_nid_dict
    neighbors = []
    while len(neighbors) < neighborhood_size:
        # prune-reattach only
        cp_tree, cp_dict = copy_tree(start_tree)
        node_ids = list(cp_dict.keys())
        prune = random.choice(node_ids)
        reattach = random.choice(node_ids)

        if prune != reattach:
            if prune_and_reattach(cp_dict[prune], cp_dict[reattach], cp_dict):
                neighbors.append((cp_tree, cp_dict))
            else:
                cp_tree = None
                cp_dict = None
    return neighbors


def hill_climbing(
    start_tree,
    start_nid_dict,
    neighborhood_size,
    max_iterations,
    alpha,
    beta,
    input_scs,
):
    current_tree = start_tree
    current_dict = start_nid_dict
    current_lh, _ = greedy_tree_likelihood(
        start_tree, start_nid_dict, input_scs, alpha, beta
    )
    tsc.logg.debug("Initial log-likelihood: %f" % current_lh)

    current_iteration = 1
    while current_iteration < max_iterations:

        if current_iteration % 10 == 0:
            tsc.logg.debug("Current iteration: %d" % current_iteration)

        neighbors = generate_neighborhood(current_tree, current_dict, neighborhood_size)
        next_eval = -np.inf
        next_sol = None

        for ng in neighbors:
            # tsc.logg.debug(ng)
            ng_lh, _ = greedy_tree_likelihood(ng[0], ng[1], input_scs, alpha, beta)

            if ng_lh > next_eval:
                next_eval = ng_lh
                next_sol = (ng[0], ng[1])

        if next_eval > current_lh:
            current_tree = next_sol[0]
            current_dict = next_sol[1]
            current_lh = next_eval
            tsc.logg.debug("Found a better solution with likelihood: %f" % current_lh)
        current_iteration += 1

    return current_tree, current_dict


def gpps_hc(input_matrix, ilp_matrix, alpha, beta, k_dollo, mut_names, ns=30, mi=100):
    imported_tree, imported_nid_dict = import_ilp_out(ilp_matrix, k_dollo, mut_names)

    hc_best_tree, hc_best_dict = hill_climbing(
        imported_tree,
        imported_nid_dict,
        neighborhood_size=ns,
        max_iterations=mi,
        alpha=alpha,
        beta=beta,
        input_scs=input_matrix,
    )

    e_mat = get_expect_matrix(hc_best_tree, hc_best_dict, input_matrix, alpha, beta)

    return e_mat
