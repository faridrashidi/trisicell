import copy
import os
import time
from decimal import *

import Bio.Phylo as bp
import numpy as np
import pandas as pd
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

import trisicell as tsc


def scistree(df_input, alpha, beta, experiment=False):
    """Solving using ScisTree.

    Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach :cite:`ScisTree`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    if not experiment:
        tsc.logg.info(f"running ScisTree with alpha={alpha}, beta={beta}")
    tmpdir = tsc.ul.tmpdirsys(suffix=".scistree")
    cells = df_input.index
    snvs = df_input.columns
    matrix_input = df_input.values
    df = df_input.transpose()

    df = df.replace(3, 0.5)
    df = df.replace(0, 1 - beta)
    df = df.replace(1, alpha)

    file1 = f"{tmpdir.name}/scistree.input"
    df.index.name = f"HAPLOID {df.shape[0]} {df.shape[1]}"
    df.to_csv(file1, sep=" ")
    with open(file1, "r") as ifile:
        data = ifile.read()
    with open(file1, "w") as ofile:
        data = data.replace('"', "")
        ofile.write(data)

    scistree = tsc.ul.get_file("trisicell.external/bin/scistree")
    cmd = (
        f"{scistree} "
        "-v "
        "-d 0 "
        "-e "
        f"-o {tmpdir.name}/scistree.gml "
        f"{tmpdir.name}/scistree.input > {tmpdir.name}/scistree.output"
    )
    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    data = []
    mut_tree = ""
    cell_tree = ""
    detail = {"cost": "\n"}
    with open(f"{tmpdir.name}/scistree.output") as infile:
        now_store = False
        for line in infile:
            line = line.strip()
            if "Mutation tree:" in line:
                mut_tree = line.split(":")[1].replace(" ", "")
                mut_tree = mut_tree.replace("#", "")
            if "Constructed single cell phylogeny:" in line:
                cell_tree = line.split(":")[1].replace(" ", "")
            if "Imputed genotypes:" in line:
                now_store = True
            if line[:4] == "Site" and now_store:
                line = "".join(line.split(":")[1])
                line = line.replace("\t", "")
                data.append([int(x) for x in line.split(" ")])
            if "current cost: " in line:
                cost = float(line.split("current cost: ")[1].split(", opt tree: ")[0])
                detail["cost"] += f"    current best cost = {cost}\n"

    data = np.array(data)
    matrix_output = data.T

    df_output = pd.DataFrame(matrix_output)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tmpdir.cleanup()

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
        # for k, v in detail.items():
        #     tsc.logg.info(f"{k}: {v}")
        return df_output
    else:
        return df_output, running_time


def rscistree(adata, mode="haploid"):
    tsc.logg.info(f"running rScisTree with mode={mode}")
    tmpdir = tsc.ul.tmpdirsys(suffix=".scistree", dirname=".")

    cells = adata.obs_names
    snvs = adata.var_names
    tsc.pp.build_scmatrix(adata)
    df_input = adata.to_df()

    V = adata.layers["mutant"]
    R = adata.layers["total"] - V
    with open(f"{tmpdir.name}/rscistree.counts", "w") as fout:
        fout.write(f"{mode.upper()} {len(snvs)} {len(cells)}\n")
        for j in range(len(snvs)):
            for i in range(len(cells)):
                fout.write(f"{R[i,j]} {V[i,j]}     ")
            fout.write(f"\n")
    cmd = f"/home/frashidi/software/temp/scistree/scprob/scprob_{mode.upper()} "
    cmd += f"{tmpdir.name}/rscistree.counts > {tmpdir.name}/rscistree.input"
    os.system(cmd)

    scistree = tsc.ul.get_file("trisicell.external/bin/scistree")
    cmd = (
        f"{scistree} "
        "-v "
        "-d 0 "
        "-e "
        f"-o {tmpdir.name}/rscistree.gml "
        f"{tmpdir.name}/rscistree.input > {tmpdir.name}/rscistree.output"
    )
    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    data = []
    mut_tree = ""
    cell_tree = ""
    detail = {"cost": "\n"}
    with open(f"{tmpdir.name}/rscistree.output") as infile:
        now_store = False
        for line in infile:
            line = line.strip()
            if "Mutation tree:" in line:
                mut_tree = line.split(":")[1].replace(" ", "")
                mut_tree = mut_tree.replace("#", "")
            if "Constructed single cell phylogeny:" in line:
                cell_tree = line.split(":")[1].replace(" ", "")
            if "Imputed genotypes:" in line:
                now_store = True
            if line[:4] == "Site" and now_store:
                line = "".join(line.split(":")[1])
                line = line.replace("\t", "")
                data.append([int(x) for x in line.split(" ")])
            if "current cost: " in line:
                cost = float(line.split("current cost: ")[1].split(", opt tree: ")[0])
                detail["cost"] += f"    current best cost = {cost}\n"

    data = np.array(data)
    matrix_output = data.T

    df_output = pd.DataFrame(matrix_output)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tmpdir.cleanup()

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)
    # for k, v in detail.items():
    #     tsc.logg.info(f"{k}: {v}")

    return df_output


def iscistree(df_input, alpha, beta, n_iters=np.inf):
    tsc.logg.info(
        f"running iScisTree with alpha={alpha}, beta={beta}, n_iters={n_iters}"
    )

    def get_initial_tree(D):
        Q = []
        for i in range(D.shape[0]):
            Q.append(list(D[i, : i + 1]))
        constructor = DistanceTreeConstructor()
        dm = DistanceMatrix(names=[f"{i}" for i in range(D.shape[0])], matrix=Q)
        tree = nj(dm)
        # tree = constructor.nj(dm)
        # tree = constructor.upgma(dm)

        node = None
        for clade in tree.find_clades():
            if clade.name == f"{D.shape[0]-1}":
                node = clade
        tree.root_with_outgroup(node)
        tree.prune(f"{D.shape[0]-1}")
        return tree

    def get_subtrees(tree):
        subtrees = []
        for clade in tree.find_clades():
            gen = np.zeros(tree.count_terminals(), dtype=int)
            cel = [
                int(clade2.name)
                for clade2 in clade.find_clades()
                if "Inner" not in clade2.name
            ]
            gen[cel] = 1
            subtrees.append(gen)
        subtrees = np.array(subtrees)
        return subtrees

    def denoise_quadratic(I, alpha, beta, subtrees):
        def column_pairs_cost(A, Ap, unit_costs):
            num = np.zeros((2, 2), dtype=np.int)
            for i in range(2):
                for j in range(2):
                    num[i, j] = np.count_nonzero(np.logical_and(A == i, Ap == j))
            return np.sum(num * unit_costs)

        unit_prob = np.array([[1 - beta, beta], [alpha, 1 - alpha]])
        unit_costs = -np.log(unit_prob)
        output = np.zeros(I.shape, dtype=int)
        total_cost = 0
        for c in range(I.shape[1]):
            costs = [
                column_pairs_cost(I[:, c], subtrees[st_ind], unit_costs)
                for st_ind in range(len(subtrees))
            ]
            ind = np.argmin(costs)
            output[:, c] = subtrees[ind]
            total_cost += costs[ind]
        return output, total_cost

    def denoise_linear(I, alpha, beta, opt_tree):
        tree = {}
        for clade in list(opt_tree.find_clades(order="level"))[::-1]:
            children = list(clade.find_clades(order="level"))
            if len(children) > 2:
                child_l = children[2]
                child_r = children[1]
                tree[clade.name] = [child_l.name, child_r.name]
            else:
                tree[clade.name] = []

        def get_cells_in_best(cells_in_best, best):
            for node in tree[best]:
                if "Inner" in node:
                    for child in tree[node]:
                        get_cells_in_best(cells_in_best, child)
                else:
                    cells_in_best.append(node)
            if "Inner" not in best:
                cells_in_best.append(best)

        output = np.zeros(I.shape, dtype=int)
        total_cost = 0
        for c in range(I.shape[1]):
            qs = {}
            best = None
            best_v = 0
            for k, v in tree.items():
                if len(v) == 0:
                    obs = I[int(k), c] == 1
                    # qs[k] = (beta**(1-obs) + (1-beta)**obs) / (alpha**obs + (1-alpha)**(1-obs))
                    # qs[k] = np.log((beta ** (1 - obs)) / (alpha ** obs))
                    p0 = (1 - obs) * (1 - beta) + obs * alpha
                    qs[k] = np.log((1 - p0) / p0)
                else:
                    qs[k] = qs[v[0]] + qs[v[1]]
                if qs[k] > best_v:
                    best = k
                    best_v = qs[k]

            cells_in_best = []
            get_cells_in_best(cells_in_best, best)
            output[list(map(int, cells_in_best)), c] = 1

            # a = (I[:, c] == 0).sum() * np.log(1 - beta)
            # b = (I[:, c] == 1).sum() * np.log(alpha)
            # best_v += a + b
            total_cost += -best_v
        return output, total_cost

    def draw_helper(x):
        if x.is_terminal():
            return f"{int(x.name)+1}"
        else:
            return None

    def get_neighbors(tree):
        """
        For a tree with n taxa, there are n - 3 internal branches.
        Thus there are 2(n - 3) NNI rearrangements for any tree
        """
        # make child to parent dict
        parents = {}
        for clade in tree.find_clades():
            if clade != tree.root:
                node_path = tree.get_path(clade)
                # cannot get the parent if the parent is root. Bug?
                if len(node_path) == 1:
                    parents[clade] = tree.root
                else:
                    parents[clade] = node_path[-2]
        neighbors = []
        root_childs = []
        for clade in tree.get_nonterminals(order="level"):
            if clade == tree.root:
                left = clade.clades[0]
                right = clade.clades[1]
                root_childs.append(left)
                root_childs.append(right)
                if not left.is_terminal() and not right.is_terminal():
                    # make changes around the left_left clade
                    # left_left = left.clades[0]
                    left_right = left.clades[1]
                    right_left = right.clades[0]
                    right_right = right.clades[1]
                    # neightbor 1 (left_left + right_right)
                    del left.clades[1]
                    del right.clades[1]
                    left.clades.append(right_right)
                    right.clades.append(left_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (left_left + right_left)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(right_left)
                    right.clades.append(right_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (left_left + left_right)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(left_right)
                    right.clades.insert(0, right_left)
            elif clade in root_childs:
                # skip root child
                continue
            else:
                # method for other clades
                # make changes around the parent clade
                left = clade.clades[0]
                right = clade.clades[1]
                parent = parents[clade]
                if clade == parent.clades[0]:
                    sister = parent.clades[1]
                    # neighbor 1 (parent + right)
                    del parent.clades[1]
                    del clade.clades[1]
                    parent.clades.append(right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(sister)
                    clade.clades.insert(0, left)
                else:
                    sister = parent.clades[0]
                    # neighbor 1 (parent + right)
                    del parent.clades[0]
                    del clade.clades[1]
                    parent.clades.insert(0, right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, sister)
                    clade.clades.insert(0, left)
        return neighbors

    def nj(distance_matrix):
        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]
            return BaseTree.Tree(root, rooted=False)
        elif len(dm) == 2:
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = BaseTree.Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]
            return BaseTree.Tree(root, rooted=False)
        while len(dm) > 2:
            # calculate nodeDist
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)

            # find minimum distance pair
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
                        min_i = i
                        min_j = j
            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = (
                dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]
            ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            if clade1.branch_length < 0:
                clade1.branch_length = 0
            if clade2.branch_length < 0:
                clade2.branch_length = 0

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (
                        dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
                    ) / 2.0

            dm.names[min_j] = "Inner" + str(inner_count)
            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return BaseTree.Tree(root, rooted=False)

    cells = list(df_input.index)
    snvs = list(df_input.columns)
    I = df_input.values

    s_time = time.time()
    Ip = np.vstack([I, np.zeros(I.shape[1])])  # add root with profile zero
    dist = tsc.ul.dist_l1_ignore_na(Ip)
    opt_tree = get_initial_tree(dist)
    # opt_subtrees = get_subtrees(opt_tree)
    # opt_O, opt_cost = denoise_quadratic(I, alpha, beta, opt_subtrees)
    opt_O, opt_cost = denoise_linear(I, alpha, beta, opt_tree)
    tsc.logg.info("current best cost =", opt_cost, time=True)

    n_iter = 1
    is_done = False
    already_seen = set()
    already_seen.add(str(opt_tree))
    while not is_done and n_iter < n_iters:
        is_done = True
        neighbors = get_neighbors(opt_tree)
        for nbr_tree in neighbors:
            if str(nbr_tree) in already_seen:
                continue
            else:
                already_seen.add(str(nbr_tree))
            # nbr_subtrees = get_subtrees(nbr_tree)
            # nbr_O, nbr_cost = denoise_quadratic(I, alpha, beta, nbr_subtrees)
            nbr_O, nbr_cost = denoise_linear(I, alpha, beta, nbr_tree)
            if nbr_cost < opt_cost:
                opt_tree = nbr_tree
                # opt_subtrees = nbr_subtrees
                opt_O = nbr_O
                opt_cost = nbr_cost
                is_done = False
                tsc.logg.info("current best cost =", nbr_cost, time=True)
        n_iter += 1
    e_time = time.time()
    running_time = e_time - s_time

    df_output = pd.DataFrame(opt_O)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
