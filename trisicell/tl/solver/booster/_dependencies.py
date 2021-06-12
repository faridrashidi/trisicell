import os

import pandas as pd
from tqdm import tqdm

DIFFERENT_LINEAGES = 0
ANCESTOR_DESCENDANT = 1
DESCENDANT_ANCESTOR = 2
SAME_NODE = 3
UNDEFINED_DEPENDENCY = 4


def read_CF_matrix_to_int_dictionary(path_CFmatrix_file):
    assert os.path.exists(path_CFmatrix_file), "ERROR. No file " + path_input_file

    input_file = open(path_CFmatrix_file)
    mut_ids = input_file.readline().strip().split()[1:]
    m = len(mut_ids)
    assert len(set(mut_ids)) == m, "ERROR. Repeating mutations present in column names."

    E = {}
    for line in input_file.readlines():
        line_columns = line.strip().split()
        assert len(line_columns) == m + 1, (
            "ERROR. Input file not in desired matrix format (note that tab or empty"
            " spaces are assumed as separators)"
        )
        cell_id = line_columns[0]
        assert cell_id not in E, "ERROR. Repeating cell id " + cell_id
        E[cell_id] = {}

        for i in range(m):
            mut_id = mut_ids[i]
            value = int(line_columns[i + 1])
            assert value == 0 or value == 1, (
                "ERROR. Conflict free matrix "
                + path_CFmatrix_file
                + " contains non-binary entries"
            )
            E[cell_id][mut_id] = value

    input_file.close()

    return E, mut_ids


def get_dependency_from_conflict_free_matrix(mut1, mut2, matrix):
    """
    This function is used to compute phylogenetic dependency between two mutations based on the conflict-free matrix.

    Arguments:
    mut1: str
        ID of the first mutation. Should be among the columns of matrix.
    mut2: str
        ID of the second mutation. Should be among the columns of matrix.
    matrix: dictionary of dictionarie
        matrix is represented as dictionary of dictionaries with elements matrix[cell_id][mut_id]

    Returns:
    -------
    This function can return one of the following five values or report error (unknown case encountered).
    ANCESTOR_DESCENDANT: if there exists at least one cell where mut1=1 and mut2=0 and at least one cell where mut1=mut2=1
    DESCENDANT_ANCESTOR: if there exists at least one cell where mut1=0 and mut2=1 and at least one cell where mut1=mut2=1
    DIFFERENT_LINEAGES: if there exists at least one cell where mut1=1 and mut2=0 and at least one cell where mut1=0 and mut2=1
    SAME_NODE: if there exists at least one cell where mut1=mut2=1 and no cells where mut1=1 and mut2=0 or mut1=0 and mut2=1.
    UNDEFINED_DEPENDENCY: (ideally) if none of the above is returned.
    """

    pairs = [[0 for i in range(2)] for j in range(2)]
    for cell_id in list(matrix.keys()):
        a = int(matrix[cell_id][mut1])
        b = int(matrix[cell_id][mut2])
        if a not in [0, 1] or b not in [0, 1]:
            continue
        else:
            pairs[a][b] += 1

    if pairs[0][1] > 0 and pairs[1][0] > 0:
        assert pairs[1][1] == 0, (
            "ERROR. There exists conflict in the input matrix for mutations "
            + mut1
            + " and "
            + mut2
        )
        return DIFFERENT_LINEAGES

    if (
        pairs[0][1] > 0
    ):  # we know that in this case pairs[1][0] = 0 (see if statement above)
        if pairs[1][1] == 0:
            return UNDEFINED_DEPENDENCY
        else:
            return DESCENDANT_ANCESTOR

    if (
        pairs[1][0] > 0
    ):  # we know that in this case pairs[0][1]  = 0 (see previous to the above if statement above)
        if pairs[1][1] == 0:
            return UNDEFINED_DEPENDENCY
        else:
            return ANCESTOR_DESCENDANT

    if pairs[1][0] == 0 and pairs[0][1] == 0:
        if pairs[1][1] > 0:
            return SAME_NODE
        else:
            return UNDEFINED_DEPENDENCY

    assert False, (
        "ERROR. Encountered case not covered in the implementation of function"
        " get_dependency."
    )


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def float_to_string(float_number, num_decimals=2):
    if is_float(float_number):
        return ("{:." + str(num_decimals) + "f}").format(float_number)
    else:
        return float_number


def prepare_dependencies(
    all_mut_ids, input_folder, output_file, max_num_submatrices, disable_tqdm
):
    assert len(set(all_mut_ids)) == len(
        all_mut_ids
    ), "ERROR. Mutation ids must be unique."

    dependencies = {}
    for mut1 in all_mut_ids:
        dependencies[mut1] = {}
        for mut2 in all_mut_ids:
            dependencies[mut1][mut2] = {}
            for dependency in [
                DIFFERENT_LINEAGES,
                ANCESTOR_DESCENDANT,
                DESCENDANT_ANCESTOR,
                SAME_NODE,
                UNDEFINED_DEPENDENCY,
            ]:
                dependencies[mut1][mut2][dependency] = 0

    CF_matrices_filenames = [
        x for x in os.listdir(input_folder) if x.endswith(".CFMatrix")
    ]
    # total = len(CF_matrices_filenames)
    # print(
    #     os.path.splitext(os.path.basename(noisy_SC_file))[0],
    #     len(CF_matrices_filenames),
    #     max_num_submatrices,
    # )

    completed = 0
    for filename in tqdm(
        CF_matrices_filenames,
        ascii=True,
        ncols=100,
        desc="DEPENDENCIES   (2/3)",
        total=max_num_submatrices,
        position=0,
        disable=disable_tqdm,
    ):
        if max_num_submatrices is not None and completed == max_num_submatrices:
            break

        # print("Completed " + str(completed) + " out of " + str(total))
        try:
            matrix, mut_ids = read_CF_matrix_to_int_dictionary(
                os.path.join(input_folder, filename)
            )
        except:
            completed -= 1

        for i in range(len(mut_ids)):
            mut1 = mut_ids[i]
            dependencies[mut1][mut1][SAME_NODE] += 1

            for j in range(i + 1, len(mut_ids)):
                mut2 = mut_ids[j]
                dependency = get_dependency_from_conflict_free_matrix(
                    mut1, mut2, matrix
                )

                if (
                    dependency == SAME_NODE
                    or dependency == DIFFERENT_LINEAGES
                    or dependency == UNDEFINED_DEPENDENCY
                ):
                    dependencies[mut1][mut2][dependency] += 1
                    dependencies[mut2][mut1][dependency] += 1
                elif dependency == ANCESTOR_DESCENDANT:
                    dependencies[mut1][mut2][ANCESTOR_DESCENDANT] += 1
                    dependencies[mut2][mut1][DESCENDANT_ANCESTOR] += 1
                elif dependency == DESCENDANT_ANCESTOR:
                    dependencies[mut1][mut2][DESCENDANT_ANCESTOR] += 1
                    dependencies[mut2][mut1][ANCESTOR_DESCENDANT] += 1
                else:
                    assert False, (
                        "ERROR. Encountered unknown dependency in file "
                        + os.path.join(input_folder, filename)
                        + " for mutations "
                        + mut1
                        + " and "
                        + mut2
                    )
        completed += 1

    output_file = open(output_file, "w")
    output_file.write(
        "\t".join(
            [
                "MUT1",
                "MUT2",
                "DIFFERENT_LINEAGES",
                "ANCESTOR_DESCENDANT",
                "DESCENDANT_ANCESTOR",
                "SAME_NODE",
                "UNDEFINED_DEPENDENCY",
                "TOTAL",
                "DOMINANT_DEPENDENCY",
                "PERCENTAGE_DOMINANT",
            ]
        )
    )
    output_file.write("\n")

    for mut1 in all_mut_ids:
        for mut2 in all_mut_ids:
            line_to_write = mut1 + "\t" + mut2

            dominant_dependency = None
            max_count = 0
            total_count = 0
            for dependency in [
                DIFFERENT_LINEAGES,
                ANCESTOR_DESCENDANT,
                DESCENDANT_ANCESTOR,
                SAME_NODE,
                UNDEFINED_DEPENDENCY,
            ]:
                current_count = dependencies[mut1][mut2][dependency]
                line_to_write += "\t" + str(current_count)

                if current_count > max_count:
                    dominant_dependency = dependency
                    max_count = current_count

                total_count += current_count

            line_to_write += "\t" + str(total_count)

            if total_count == 0:
                line_to_write += "\t" + "UNDEFINED"
                line_to_write += "\t" + "UNDEFINED"
                if mut1 < mut2:
                    # print("WARNING. For file " + noisy_SC_file + " mutations " + mut1 + " and " + mut2 + " were never subsampled together")
                    pass
            else:
                line_to_write += "\t" + str(dominant_dependency)
                line_to_write += "\t" + float_to_string(
                    (100.0 * max_count) / total_count
                )

            output_file.write(line_to_write + "\n")

    output_file.close()
