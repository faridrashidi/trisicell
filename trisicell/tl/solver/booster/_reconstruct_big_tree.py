import argparse
import csv
import math
import os

import numpy as np
from tqdm import tqdm

import trisicell as tsc

DIFFERENT_LINEAGES = 0
ANCESTOR_DESCENDANT = 1
DESCENDANT_ANCESTOR = 2
SAME_NODE = 3
UNDEFINED_DEPENDENCY = 4
ROOT_NODE_ID = "-1"  # this is ID of the root node in tree representation
ROOT = ROOT_NODE_ID


def read_matrix_from_file_into_hash(path_matrix_file):
    """
    This function reads tab or space delimited file into 2D-dictionary.

    TO DO
    -----
        Add full documentation of the function (not needed at the moment)
    """

    assert os.path.exists(path_matrix_file), (
        "ERROR. There does not exist file " + path_matrix_file
    )

    input_file = open(path_matrix_file)
    input_lines = input_file.readlines()
    input_file.close()

    assert len(input_lines) > 0, "ERROR. Input file " + path_matrix_file + " is empty."

    header_entries = input_lines[0].strip().split()
    column_ids = header_entries[1 : len(header_entries)]
    num_columns = len(column_ids)
    assert num_columns != 0, (
        "ERROR. Number of columns in " + path_matrix_file + " equals zero. Exiting!!!"
    )

    input_lines_without_header = input_lines[1 : len(input_lines)]
    D = {}
    for line in input_lines_without_header:
        line_columns = line.strip().split()
        row_id = line_columns[0].strip()
        assert row_id not in D.keys(), "ERROR. " + row_id + " is already in keys."
        D[row_id] = {}
        for i in range(num_columns):
            assert column_ids[i] not in D[row_id].keys(), (
                "ERROR. " + column_ids[i] + " is already in keys."
            )
            D[row_id][column_ids[i]] = line_columns[i + 1]
    return D


def hash_entries_to_integers(D):
    """
    This function simply converts all entries of a given dictionary to integers.

    Arguments:
    ---------
    D: dict
        Dictionary of dictionaries (2D-hash)

    Errors:
    ------
        not covered in the current implementation
    """
    for row in D:
        for column in D[row]:
            D[row][column] = int(D[row][column])
    return D


def write_dictionary_of_dictionaries_to_file(D, path_output_file):
    """
    The function below can be used to write dictionary of dictionaries D to a file.
    D is in most applications single-cell matrix, e.g. D[cell_id][mut_id] = 1

    Arguments:
    ---------
    D: dictionary of dictionaries
    path_output_file: String (if file exists it will be overwritten)

    """

    row_ids = list(D.keys())
    column_ids = list(D[row_ids[0]].keys())
    assert (
        len(row_ids) > 0 and len(column_ids) > 0
    ), "ERROR. Number of columns or number of rows equal to zero."

    for row_id in row_ids:
        assert set(D[row_id].keys()) == set(column_ids), (
            "ERROR. Some rows have different column ids so this hash is not in a matrix"
            " form,"
        )

    output_file = open(path_output_file, "w")

    output_file.write("row_x_col" + "\t" + "\t".join(column_ids) + "\n")

    for row_id in row_ids:
        output_file.write(
            row_id
            + "\t"
            + "\t".join([str(D[row_id][column_id]) for column_id in column_ids])
            + "\n"
        )

    output_file.close()


def get_row_ids_from_2D_hash(D):
    """
    This function returns IDs of rows of 2-dimensional matrix stored via dictionary of dictionaries.

    Arguments:
    ---------
    D: dictionary of dictionaries

    Returns:
    -------
    Array of strings representing list of keys of D.
    """

    return list(D.keys())


def get_column_ids_from_2D_hash(D):
    r"""
    This function returns IDs of columns of 2-dimensional matrix stored via dictionary of dictionaries.

    Arguments:
    ---------
    D: dictionary of dictionaries
        It is assumed that this dictionary is such that D[x].keys() is the same for all x \in D.keys().

    Returns:
    -------
    Array of strings representing list of keys of D.
    """

    row_ids = get_row_ids_from_2D_hash(D)
    assert len(row_ids) > 0, "ERROR. Dictionary empty."
    row_0_id = row_ids[0]
    return list(D[row_0_id].keys())


def get_CF_matrix_from_parent_vector(parent, D, alpha, beta):
    """
    Documentation to be added.
    """

    cell_ids = list(D.keys())
    mut_ids = list(D[cell_ids[0]].keys())

    children = {}
    children[ROOT] = []
    for mut_id in mut_ids:
        children[mut_id] = []

    for child_id, parent_id in parent.items():
        if child_id != ROOT:
            children[parent_id].append(child_id)

    E = {}
    for cell_id in cell_ids:
        E[cell_id] = {}
        for mut_id in mut_ids:
            E[cell_id][mut_id] = None

    for cell_id in cell_ids:
        score = {}
        score[ROOT] = 0
        for mut_id in mut_ids:
            observed = int(D[cell_id][mut_id])
            if observed == 0:
                score[ROOT] += math.log(1 - alpha)
            if observed == 1:
                score[ROOT] += math.log(alpha)
        best_score = score[ROOT]
        best_mut = ROOT

        muts_to_visit = [child_id for child_id in children[ROOT]]
        while len(muts_to_visit) > 0:
            mut_id = muts_to_visit.pop(0)
            parent_id = parent[mut_id]

            score[mut_id] = score[
                parent_id
            ]  # this is only temporary. see changes below
            observed = int(D[cell_id][mut_id])
            if observed == 0:
                score[mut_id] -= math.log(1 - alpha)
                score[mut_id] += math.log(beta)
            if observed == 1:
                score[mut_id] -= math.log(alpha)
                score[mut_id] += math.log(1 - beta)

            if score[mut_id] > best_score:
                best_score = score[mut_id]
                best_mut = mut_id

            for child_id in children[mut_id]:
                muts_to_visit.append(child_id)

        muts_present_in_true_genotype = []
        current_mut = best_mut
        while current_mut != ROOT:
            muts_present_in_true_genotype.append(current_mut)
            current_mut = parent[current_mut]

        for mut_id in mut_ids:
            if mut_id in muts_present_in_true_genotype:
                E[cell_id][mut_id] = 1
            else:
                E[cell_id][mut_id] = 0

    zero_one_flips = 0
    one_zero_flips = 0

    for cell_id in cell_ids:
        for mut_id in mut_ids:
            observed = int(D[cell_id][mut_id])
            true = int(E[cell_id][mut_id])
            if observed == 1 and true == 0:
                one_zero_flips += 1
            if observed == 0 and true == 1:
                zero_one_flips += 1

    # print("0_1_flips:    " + str(zero_one_flips))
    # print("1_0_flips:    " + str(one_zero_flips))
    return E


def compute_weights_from_dependencies_file(path_dependencies_file):
    """
    This function takes as input path to dependencies file and returns weights.

    Arguments:
    ---------
    path_dependencies_file: str
        Path to file which stores dependencies counts obtained after processing folder with .CFMatrix files obtained on subsampled instances.
        Expected header of this file can be found below in expected_header_line variable and it explains columns present in dependencies file.

    Returns:
    -------
    3-D dictionary weights accessed as follows: weights[mut1][mut2][dependency]
    Here, dependency can be any of ANCESTOR_DESCENDANT, DESCENDANT_ANCESTOR or DIFFERENT_LINEAGES.
    Dictionary is not ragged (i.e., both [mut1][mut2] and [mut2][mut1] are populated).
    """

    assert os.path.exists(path_dependencies_file), (
        "ERROR. There does not exist file " + path_dependencies_file
    )

    dependencies_file = open(path_dependencies_file)
    expected_header_line = "MUT1	MUT2	DIFFERENT_LINEAGES	ANCESTOR_DESCENDANT	DESCENDANT_ANCESTOR	SAME_NODE	"
    expected_header_line += (
        "UNDEFINED_DEPENDENCY	TOTAL	DOMINANT_DEPENDENCY	PERCENTAGE_DOMINANT"
    )
    observed_header_line = dependencies_file.readline().rstrip()
    assert observed_header_line == expected_header_line, (
        "ERROR for dependencies file "
        + path_dependencies_file
        + ". Header line different from expected."
    )
    dependencies_file.close()

    dependencies_file_rows = list(
        csv.DictReader(open(path_dependencies_file), delimiter="\t")
    )

    mut_ids = list({row["MUT1"] for row in dependencies_file_rows})
    m = len(mut_ids)

    weights = {}
    weights[ROOT] = {}
    for mut_id in mut_ids:
        weights[mut_id] = {}
        weights[ROOT][mut_id] = {ANCESTOR_DESCENDANT: 0, DIFFERENT_LINEAGES: 0}
        weights[mut_id][ROOT] = {ANCESTOR_DESCENDANT: 0, DIFFERENT_LINEAGES: 0}

    for row in dependencies_file_rows:
        mut1 = row["MUT1"]
        mut2 = row["MUT2"]

        weights[mut1][mut2] = {}

        total = (
            int(row["ANCESTOR_DESCENDANT"])
            + int(row["DESCENDANT_ANCESTOR"])
            + int(row["DIFFERENT_LINEAGES"])
            + int(row["SAME_NODE"])
        )

        if total == 0:
            weights[mut1][mut2][DIFFERENT_LINEAGES] = 0
            weights[mut1][mut2][ANCESTOR_DESCENDANT] = 0
            weights[mut1][mut2][DESCENDANT_ANCESTOR] = 0
        else:
            weights[mut1][mut2][DIFFERENT_LINEAGES] = (
                float(row["DIFFERENT_LINEAGES"])
            ) / total
            weights[mut1][mut2][ANCESTOR_DESCENDANT] = (
                float(row["ANCESTOR_DESCENDANT"]) + 0.5 * float(row["SAME_NODE"])
            ) / total
            weights[mut1][mut2][DESCENDANT_ANCESTOR] = (
                float(row["DESCENDANT_ANCESTOR"]) + 0.5 * float(row["SAME_NODE"])
            ) / total

        if row["DOMINANT_DEPENDENCY"] == "UNDEFINED":
            weights[mut1][mut2]["DOMINANT_DEPENDENCY"] = "UNDEFINED"
        else:
            weights[mut1][mut2]["DOMINANT_DEPENDENCY"] = int(row["DOMINANT_DEPENDENCY"])

        weights[mut1][mut2]["TOTAL"] = int(row["TOTAL"])

    return weights


# this function is used to give an order of mutations in which they will be added to the big tree
def get_ordered_list_of_mutations_from_input_matrix(path_noisy_SC_matrix):
    D = hash_entries_to_integers(read_matrix_from_file_into_hash(path_noisy_SC_matrix))
    cell_ids = list(D.keys())
    mut_ids = list(D[cell_ids[0]].keys())

    negative_num_ones = (
        {}
    )  # = (-1)*(number of ones in a column). Consequently, mutations having more ones will be added to tree first.
    for mut_id in mut_ids:
        negative_num_ones[mut_id] = (-1) * len(
            [D[cell_id][mut_id] for cell_id in cell_ids if D[cell_id][mut_id] == 1]
        )

    mut_ids.sort(key=lambda mut_id: negative_num_ones[mut_id])

    return mut_ids


def get_next_mutation(weights, existing_mutations):
    mut_ids = list(weights.keys())
    assert len(existing_mutations) < len(
        mut_ids
    ), "ERROR in calling function get_next_mutation. All mutations have been used."

    candidates_for_addition = [x for x in mut_ids if x not in existing_mutations]

    if len(candidates_for_addition) == 1:
        return candidates_for_addition[0]

    else:
        best_score = -1
        best_mutation = None
        for pivot in candidates_for_addition:
            current_score = np.sum(
                [
                    (
                        weights[pivot][mut_id][ANCESTOR_DESCENDANT]
                        - weights[pivot][mut_id][DESCENDANT_ANCESTOR]
                    )
                    for mut_id in candidates_for_addition
                    if mut_id != pivot
                ]
            )
            if current_score > best_score:
                best_score = current_score
                best_mutation = pivot

        return best_mutation


def get_CF_matrix_from_parent_vector(parent, D, alpha, beta):
    cell_ids = get_row_ids_from_2D_hash(D)
    mut_ids = get_column_ids_from_2D_hash(D)

    children = {}
    children[ROOT] = []
    for mut_id in mut_ids:
        children[mut_id] = []

    for child_id, parent_id in parent.items():
        if child_id != ROOT:
            children[parent_id].append(child_id)

    E = {}
    for cell_id in cell_ids:
        E[cell_id] = {}
        for mut_id in mut_ids:
            E[cell_id][mut_id] = None

    for cell_id in cell_ids:
        score = {}
        score[ROOT] = 0
        for mut_id in mut_ids:
            observed = D[cell_id][mut_id]
            if observed == 0:
                score[ROOT] += math.log(1 - alpha)
            if observed == 1:
                score[ROOT] += math.log(alpha)
        best_score = score[ROOT]
        best_mut = ROOT

        muts_to_visit = [child_id for child_id in children[ROOT]]
        while len(muts_to_visit) > 0:
            mut_id = muts_to_visit.pop(0)
            parent_id = parent[mut_id]

            score[mut_id] = score[
                parent_id
            ]  # this is only temporary. see changes below
            observed = D[cell_id][mut_id]
            if observed == 0:
                score[mut_id] -= math.log(1 - alpha)
                score[mut_id] += math.log(beta)
            if observed == 1:
                score[mut_id] -= math.log(alpha)
                score[mut_id] += math.log(1 - beta)

            if score[mut_id] > best_score:
                best_score = score[mut_id]
                best_mut = mut_id

            for child_id in children[mut_id]:
                muts_to_visit.append(child_id)

        muts_present_in_true_genotype = []
        current_mut = best_mut
        while current_mut != ROOT:
            muts_present_in_true_genotype.append(current_mut)
            current_mut = parent[current_mut]

        for mut_id in mut_ids:
            if mut_id in muts_present_in_true_genotype:
                E[cell_id][mut_id] = 1
            else:
                E[cell_id][mut_id] = 0

    return E


def reconstruct_big_tree(
    path_dependencies_file,
    path_noisy_SC_matrix,
    alpha,
    beta,
    output_files_prefix,
    disable_tqdm,
):

    weights = compute_weights_from_dependencies_file(path_dependencies_file)
    D = hash_entries_to_integers(read_matrix_from_file_into_hash(path_noisy_SC_matrix))

    ordered_mut_ids = get_ordered_list_of_mutations_from_input_matrix(
        path_noisy_SC_matrix
    )  # mutations in order that they should be added to the tree

    A = {}  # A[mut1][mut2] = 1 if and only if mut1 is ancestor of mut2
    A[ROOT] = {}
    A[ROOT][ROOT] = 0
    for mut1 in ordered_mut_ids:
        A[mut1] = {}
        A[ROOT][mut1] = 1
        A[mut1][ROOT] = 0
        for mut2 in ordered_mut_ids:
            if mut2 != mut1:
                A[mut1][mut2] = None

    parent = {}  # parent vector (keys are mutation IDs)
    parent[ROOT] = "NA"
    existing_mutations = [ROOT]  # mutations added so far to the tree
    big_tree_score = 0

    m = len(ordered_mut_ids)

    for mutation_to_add in tqdm(
        ordered_mut_ids,
        ascii=True,
        ncols=100,
        desc="RECONSTRUCTION (3/3)",
        position=0,
        disable=disable_tqdm,
    ):
        # print(
        #     "Number of mutations added: "
        #     + str(len(existing_mutations))
        #     + " out of "
        #     + str(m)
        # )

        best_choice_move = None
        best_choice_mutation = None
        best_choice_tree_score = -1

        for pivot_mut in existing_mutations:
            current_move = "child"  # add new mutation as child to the pivot mutation
            current_tree_score = big_tree_score
            for mut_id in existing_mutations:
                if (A[mut_id][pivot_mut] == 1) or (mut_id == pivot_mut):
                    current_tree_score += weights[mut_id][mutation_to_add][
                        ANCESTOR_DESCENDANT
                    ]
                else:
                    current_tree_score += weights[mut_id][mutation_to_add][
                        DIFFERENT_LINEAGES
                    ]

            if current_tree_score > best_choice_tree_score:
                best_choice_tree_score = current_tree_score
                best_choice_mutation = pivot_mut
                best_choice_move = current_move

            current_move = (  # edge (parent, pivot_mut) is broken into (parent, pivot_mut) + (pivot_mut, mutation_to_add)
                "break edge lower"
            )
            current_tree_score = big_tree_score
            for mut_id in existing_mutations:
                if (A[mut_id][pivot_mut] == 1) or (mut_id == pivot_mut):
                    current_tree_score += weights[mut_id][mutation_to_add][
                        ANCESTOR_DESCENDANT
                    ]
                elif A[pivot_mut][mut_id] == 1:
                    current_tree_score += weights[mutation_to_add][mut_id][
                        ANCESTOR_DESCENDANT
                    ]
                else:
                    current_tree_score += weights[mut_id][mutation_to_add][
                        DIFFERENT_LINEAGES
                    ]

            if current_tree_score > best_choice_tree_score:
                best_choice_tree_score = current_tree_score
                best_choice_mutation = pivot_mut
                best_choice_move = current_move

            current_move = (  # edge (parent, pivot_mut) is broken into (parent, mutation_to_add) + (mutation_to_add, pivot_mut)
                "break edge middle"
            )
            if pivot_mut == ROOT:
                continue
            # we just use previous current score with swapping pivot and mutation_to_add
            current_tree_score -= weights[pivot_mut][mutation_to_add][
                ANCESTOR_DESCENDANT
            ]
            current_tree_score += weights[mutation_to_add][pivot_mut][
                ANCESTOR_DESCENDANT
            ]
            if current_tree_score > best_choice_tree_score:
                best_choice_tree_score = current_tree_score
                best_choice_mutation = pivot_mut
                best_choice_move = current_move

        big_tree_score = best_choice_tree_score
        if best_choice_move == "child":
            parent[mutation_to_add] = best_choice_mutation
            for mut_id in existing_mutations:
                if (
                    A[mut_id][best_choice_mutation] == 1
                    or mut_id == best_choice_mutation
                ):
                    A[mut_id][mutation_to_add] = 1
                    A[mutation_to_add][mut_id] = 0
                else:
                    A[mut_id][mutation_to_add] = 0
                    A[mutation_to_add][mut_id] = 0

        if best_choice_move == "break edge lower":
            for mut_id in existing_mutations:
                if parent[mut_id] == best_choice_mutation:
                    parent[mut_id] = mutation_to_add
            parent[mutation_to_add] = best_choice_mutation

        if best_choice_move == "break edge middle":
            parent[mutation_to_add] = parent[best_choice_mutation]
            parent[best_choice_mutation] = mutation_to_add

        if best_choice_move in ["break edge lower", "break edge middle"]:
            for mut_id in existing_mutations:
                if mut_id == best_choice_mutation:
                    if best_choice_move == "break edge lower":
                        A[mut_id][mutation_to_add] = 1
                        A[mutation_to_add][mut_id] = 0
                    else:
                        A[mut_id][mutation_to_add] = 0
                        A[mutation_to_add][mut_id] = 1
                elif A[best_choice_mutation][mut_id] == 1:
                    A[mutation_to_add][mut_id] = 1
                    A[mut_id][mutation_to_add] = 0
                elif A[mut_id][best_choice_mutation] == 1:
                    A[mut_id][mutation_to_add] = 1
                    A[mutation_to_add][mut_id] = 0
                else:
                    A[mutation_to_add][mut_id] = 0
                    A[mut_id][mutation_to_add] = 0

        existing_mutations.append(mutation_to_add)
        A[mutation_to_add][mutation_to_add] = 0

    parent_vector_file = open(output_files_prefix + ".dnc.tree.parent", "w")
    for mut_id, parent_id in parent.items():
        parent_vector_file.write(mut_id + "\t" + parent_id + "\n")
    parent_vector_file.close()

    E = get_CF_matrix_from_parent_vector(parent, D, alpha, beta)
    write_dictionary_of_dictionaries_to_file(E, output_files_prefix + ".dnc.CFMatrix")

    return big_tree_score
