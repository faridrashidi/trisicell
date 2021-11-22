import copy
import time

import numpy as np
import pandas as pd
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceMatrix

import trisicell as tsc
from trisicell.external._scistree import run_scistree
from trisicell.external._scprob import run_scprob

# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


def scistree(df_input, alpha, beta, n_threads=1, experiment=False):
    """Solving using ScisTree.

    Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach
    :cite:`ScisTree`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    n_threads : :obj:`int`
        Number of threads.
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    if not experiment:
        tsc.logg.info(
            f"running ScisTree with alpha={alpha}, beta={beta}, n_threads={n_threads}"
        )
    tmpdir = tsc.ul.tmpdirsys(suffix=".scistree")
    cells = df_input.index
    snvs = df_input.columns
    df = df_input.transpose()

    df = df.replace(3, 0.5)
    df = df.replace(0, 1 - beta)
    df = df.replace(1, alpha)

    file1 = f"{tmpdir.name}/scistree.input"
    df.index.name = f"HAPLOID {df.shape[0]} {df.shape[1]}"
    df.to_csv(file1, sep=" ")
    with open(file1) as ifile:
        data = ifile.read()
    with open(file1, "w") as ofile:
        data = data.replace('"', "")
        ofile.write(data)

    cmd = [
        "scistree",
        "-v",
        "-d",
        "0",
        "-e",
        "-k",
        f"{n_threads}",
        "-o",
        f"{tmpdir.name}/scistree.gml",
        f"{tmpdir.name}/scistree.input",
    ]
    s_time = time.time()
    run_scistree(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    data = []
    detail = {"cost": "\n"}
    with open(f"{tmpdir.name}/scistree.output") as infile:
        now_store = False
        for line in infile:
            line = line.strip()
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
        return df_output
    else:
        return df_output, running_time


def rscistree(adata, alpha=0, beta=0, n_threads=1, mode="haploid"):
    """Solving using read-count ScisTree.

    Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach
    :cite:`ScisTree`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    n_threads : :obj:`int`
        Number of threads.
    mode : :obj:`str`
        Mode of calculating the probability from read-count.
        In {'haploid', 'ternary'}, by default haploid
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    tsc.logg.info(f"running rScisTree with n_threads={n_threads}, mode={mode}")
    tmpdir = tsc.ul.tmpdirsys(suffix=".rscistree", dirname=".")

    cells = adata.obs_names
    snvs = adata.var_names
    df_input = adata.to_df()

    V = adata.layers["mutant"]
    R = adata.layers["total"] - V
    with open(f"{tmpdir.name}/rscistree.counts", "w") as fout:
        fout.write(f"{mode.upper()} {len(snvs)} {len(cells)}\n")
        for j in range(len(snvs)):
            for i in range(len(cells)):
                fout.write(f"{R[i,j]} {V[i,j]}     ")
            fout.write("\n")

    cmd = [
        "scprob",
        f"{tmpdir.name}/rscistree.counts",
    ]
    if mode.lower() == "haploid":
        cmd += ["0"]
    elif mode.lower() == "ternary":
        cmd += ["1"]
    else:
        tsc.logg.error("Wrong mode!")
    run_scprob(cmd)

    cmd = [
        "scistree",
        "-v",
        "-d",
        "0",
        "-e",
        "-k",
        f"{n_threads}",
        "-o",
        f"{tmpdir.name}/rscistree.gml",
        f"{tmpdir.name}/rscistree.input",
    ]
    s_time = time.time()
    run_scistree(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    data = []
    detail = {"cost": "\n"}
    with open(f"{tmpdir.name}/rscistree.output") as infile:
        now_store = False
        for line in infile:
            line = line.strip()
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
    return df_output


def iscistree(df_input, alpha, beta, n_iters=np.inf):
    """Solving using my own implementation of ScisTree.

    Accurate and efficient cell lineage tree inference from noisy
    single cell data: the maximum likelihood perfect phylogeny approach
    :cite:`ScisTree`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    n_iters : :obj:`int`
        Number of iterations to search for the neighboring trees, by default inf.
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    tsc.logg.info(
        f"running iScisTree with alpha={alpha}, beta={beta}, n_iters={n_iters}"
    )

    def get_initial_tree(D):
        Q = []
        for i in range(D.shape[0]):
            Q.append(list(D[i, : i + 1]))
        dm = DistanceMatrix(names=[f"{i}" for i in range(D.shape[0])], matrix=Q)
        tree = nj(dm)

        node = None
        for clade in tree.find_clades():
            if clade.name == f"{D.shape[0]-1}":
                node = clade
        tree.root_with_outgroup(node)
        tree.prune(f"{D.shape[0]-1}")
        return tree

    def denoise_linear(I_mtr, alpha, beta, opt_tree):
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

        output = np.zeros(I_mtr.shape, dtype=int)
        total_cost = 0
        for c in range(I_mtr.shape[1]):
            qs = {}
            best = None
            best_v = 0
            for k, v in tree.items():
                if len(v) == 0:
                    obs = I_mtr[int(k), c] == 1
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

            total_cost += -best_v
        return output, total_cost

    def get_neighbors(tree):
        """
        Return neighbors.

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
    I_mtr = df_input.values

    s_time = time.time()
    Ip = np.vstack([I_mtr, np.zeros(I_mtr.shape[1])])  # add root with profile zero
    tsc.logg.debug("now calculating distance!", time=True)
    dist = tsc.ul.dist_l1_ignore_na(Ip)
    tsc.logg.debug("distance is done!", time=True)
    opt_tree = get_initial_tree(dist)
    # opt_subtrees = get_subtrees(opt_tree)
    # opt_O, opt_cost = denoise_quadratic(I, alpha, beta, opt_subtrees)
    tsc.logg.debug("init tree!", time=True)
    opt_O, opt_cost = denoise_linear(I_mtr, alpha, beta, opt_tree)
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
            nbr_O, nbr_cost = denoise_linear(I_mtr, alpha, beta, nbr_tree)
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
