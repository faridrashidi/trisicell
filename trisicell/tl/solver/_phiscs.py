import math
import time
from itertools import combinations

import numpy as np
import pandas as pd
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
from scipy.stats import binom

import trisicell as tsc
from trisicell.external._betabinom import pmf_BetaBinomial


def phiscsb(df_input, alpha, beta, experiment=False):
    """Solving using PhISCS-B (only in single-cell mode, no bulk).

    a combinatorial approach for subperfect tumor phylogeny reconstruction
    via integrative use of single-cell and bulk sequencing data :cite:`PhISCS`.

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
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    if not experiment:
        tsc.logg.info(f"running PhISCS-B with alpha={alpha}, beta={beta}")
    cells = list(df_input.index)
    snvs = list(df_input.columns)
    df_input = df_input.replace("?", 3)
    df_input = df_input.astype(int)
    I_mtr = df_input.values

    rc2 = RC2(WCNF())

    num_cells = len(cells)
    num_mutations = len(snvs)

    Y = np.empty((num_cells, num_mutations), dtype=np.int64)
    numVarY = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            numVarY += 1
            Y[i, j] = numVarY

    B = np.empty((num_mutations, num_mutations, 2, 2), dtype=np.int64)
    numVarB = 0
    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            for i, j in [(0, 1), (1, 0), (1, 1)]:
                numVarB += 1
                B[p, q, i, j] = numVarY + numVarB

    Z = np.empty((num_cells, num_mutations), dtype=np.int64)
    numVarZ = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            if I_mtr[i, j] == 0:
                numVarZ += 1
                Z[i, j] = numVarY + numVarB + numVarZ

    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            rc2.add_clause([-B[p, q, 0, 1], -B[p, q, 1, 0], -B[p, q, 1, 1]])
            for i in range(num_cells):
                rc2.add_clause([-Y[i, p], -Y[i, q], B[p, q, 1, 1]])
                rc2.add_clause([Y[i, p], -Y[i, q], B[p, q, 0, 1]])
                rc2.add_clause([-Y[i, p], Y[i, q], B[p, q, 1, 0]])

    for i in range(num_cells):
        for j in range(num_mutations):
            # 0->1
            if alpha == 0:
                if I_mtr[i, j] == 0:
                    rc2.add_clause([-Y[i, j]], weight=1)
                if I_mtr[i, j] == 1:
                    rc2.add_clause([Y[i, j]])

            # 0->1 and 1->0
            if alpha > 0:
                if I_mtr[i, j] == 0:
                    rc2.add_clause([Y[i, j], Z[i, j]])
                    rc2.add_clause([-Y[i, j], -Z[i, j]])
                    rc2.add_clause([Z[i, j]], weight=math.log((1 - alpha) / beta))
                if I_mtr[i, j] == 1:
                    rc2.add_clause([Y[i, j]], weight=math.log((1 - beta) / alpha))

    s_time = time.time()
    variables = rc2.compute()
    e_time = time.time()
    running_time = e_time - s_time

    sol_Y = np.empty((num_cells, num_mutations), dtype=np.int8)
    numVar = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            sol_Y[i, j] = variables[numVar] > 0
            numVar += 1

    df_output = pd.DataFrame(sol_Y)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    if not experiment:
        tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output


def phiscsi(df_input, alpha, beta, time_limit=86400, n_threads=1):
    """Solving using PhISCS-I (only in single-cell mode, no bulk).

    a combinatorial approach for subperfect tumor phylogeny reconstruction
    via integrative use of single-cell and bulk sequencing data :cite:`PhISCS`.

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
    time_limit : :obj:`int`, optional
        Time limit of the Gurobi running in seconds, by default 86400 (one day)
    n_threads : :obj:`int`, optional
        Number of threads for Gurobi solver, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    gp, gp_is_not_imported = tsc.ul.import_gurobi()
    if gp_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    tsc.logg.info(
        f"running PhISCS-I with alpha={alpha}, beta={beta}, time_limit={time_limit}, "
        f"n_threads={n_threads}"
    )
    cells = list(df_input.index)
    snvs = list(df_input.columns)
    df_input = df_input.replace("?", 3)
    df_input = df_input.astype(int)
    I_mtr = df_input.values

    model = gp.Model("ILP")
    model.Params.OutputFlag = 0
    model.Params.LogFile = ""
    model.Params.Threads = n_threads
    model.Params.TimeLimit = time_limit

    num_cells = len(cells)
    num_mutations = len(snvs)
    Y = {}
    B = {}
    for c in range(num_cells):
        for m in range(num_mutations):
            if alpha == 0:
                # 0->1
                if I_mtr[c, m] == 0:
                    Y[c, m] = model.addVar(
                        vtype=gp.GRB.BINARY, obj=1, name=f"Y({c},{m})"
                    )
                elif I_mtr[c, m] == 1:
                    Y[c, m] = 1
                else:
                    Y[c, m] = model.addVar(
                        vtype=gp.GRB.BINARY, obj=0, name=f"Y({c},{m})"
                    )
            else:
                # 0->1 & 1->0
                Y[c, m] = model.addVar(vtype=gp.GRB.BINARY, name=f"Y({c},{m})")
    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            B[p, q, 1, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,1]"
            )
            B[p, q, 1, 0] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,0]"
            )
            B[p, q, 0, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},0,1]"
            )
    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            model.addConstr(B[p, q, 0, 1] + B[p, q, 1, 0] + B[p, q, 1, 1] <= 2)
            for i in range(num_cells):
                model.addConstr(Y[i, p] + Y[i, q] - B[p, q, 1, 1] <= 1)
                model.addConstr(-Y[i, p] + Y[i, q] - B[p, q, 0, 1] <= 0)
                model.addConstr(Y[i, p] - Y[i, q] - B[p, q, 1, 0] <= 0)

    if alpha == 0:
        # 0->1
        model.Params.ModelSense = gp.GRB.MINIMIZE
    else:
        # 0->1 & 1->0
        objective = 0
        for j in range(num_mutations):
            for i in range(num_cells):
                if I_mtr[i, j] == 0:
                    objective += (
                        np.log(1 - alpha) + np.log(beta / (1 - alpha)) * Y[i, j]
                    )
                if I_mtr[i, j] == 1:
                    objective += np.log(alpha) + np.log((1 - beta) / alpha) * Y[i, j]
        model.setObjective(objective, gp.GRB.MAXIMIZE)

    s_time = time.time()
    model.optimize()
    e_time = time.time()
    running_time = e_time - s_time

    sol_Y = np.zeros((num_cells, num_mutations), dtype=np.int8)
    for i in range(num_cells):
        for j in range(num_mutations):
            if alpha == 0:
                if I_mtr[i, j] != 1:
                    sol_Y[i, j] = Y[i, j].X > 0.5
                else:
                    sol_Y[i, j] = 1
            else:
                sol_Y[i, j] = Y[i, j].X > 0.5

    df_output = pd.DataFrame(sol_Y)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output


def phiscsb_bulk(
    df_input,
    alpha,
    beta,
    kmax=0,
    vaf_info=None,
    delta=0.2,
):
    """Solving using PhISCS-B (in single-cell mode with bulk and mutation elimination).

    a combinatorial approach for subperfect tumor phylogeny reconstruction
    via integrative use of single-cell and bulk sequencing data :cite:`PhISCS`.

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
    kmax : :obj:`int`, optional
        Max number of mutations to be eliminated, by default 0
    vaf_info : :class:`pandas.DataFrame`, optional
        Information about the variant allele frequency in bulk data
        The size is n_SNVs x n_samples, by default None
    delta : :obj:`float`, optional
        Delta parameter accounting for VAF variance, by default 0.2

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    # TODO: implement

    tsc.logg.info(
        f"running PhISCS-B-orig with alpha={alpha}, beta={beta}, kmax={kmax}, "
        f"vaf_info={vaf_info}, delta={delta}"
    )

    cells = list(df_input.index)
    snvs = list(df_input.columns)
    df_input = df_input.replace("?", 3)
    df_input = df_input.astype(int)
    I_mtr = df_input.values

    rc2 = RC2(WCNF())

    num_cells = len(cells)
    num_mutations = len(snvs)

    Y = np.empty((num_cells, num_mutations), dtype=np.int64)
    numVarY = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            numVarY += 1
            Y[i, j] = numVarY

    B = np.empty((num_mutations, num_mutations, 2, 2), dtype=np.int64)
    numVarB = 0
    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            for i, j in [(0, 1), (1, 0), (1, 1)]:
                numVarB += 1
                B[p, q, i, j] = numVarY + numVarB

    Z = np.empty((num_cells, num_mutations), dtype=np.int64)
    numVarZ = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            if I_mtr[i, j] == 0:
                numVarZ += 1
                Z[i, j] = numVarY + numVarB + numVarZ

    if kmax > 0:
        K = np.empty(num_mutations + 1, dtype=np.int64)
        numVarK = 0
        for j in range(num_mutations + 1):
            numVarK += 1
            K[j] = numVarY + numVarB + numVarZ + numVarK

    if vaf_info is not None:
        A = np.empty((num_cells + 1, num_mutations + 1), dtype=np.int64)
        numVarA = 0
        for i in range(num_cells + 1):
            for j in range(num_mutations + 1):
                numVarA += 1
                A[i, j] = numVarY + numVarB + numVarZ + numVarK + numVarA

    for p in range(num_mutations):
        for q in range(p + 1, num_mutations):
            if kmax > 0:
                rc2.add_clause(
                    [K[p], K[q], -B[p, q, 0, 1], -B[p, q, 1, 0], -B[p, q, 1, 1]]
                )
            else:
                rc2.add_clause([-B[p, q, 0, 1], -B[p, q, 1, 0], -B[p, q, 1, 1]])
            for i in range(num_cells):
                rc2.add_clause([-Y[i, p], -Y[i, q], B[p, q, 1, 1]])
                rc2.add_clause([Y[i, p], -Y[i, q], B[p, q, 0, 1]])
                rc2.add_clause([-Y[i, p], Y[i, q], B[p, q, 1, 0]])

    for i in range(num_cells):
        for j in range(num_mutations):
            # 0->1
            if alpha == 0:
                if I_mtr[i, j] == 0:
                    rc2.add_clause([-Y[i, j]], weight=1)
                if I_mtr[i, j] == 1:
                    rc2.add_clause([Y[i, j]])

            # 0->1 and 1->0
            if alpha > 0:
                if I_mtr[i, j] == 0:
                    rc2.add_clause([Y[i, j], Z[i, j]])
                    rc2.add_clause([-Y[i, j], -Z[i, j]])
                    rc2.add_clause([Z[i, j]], weight=math.log((1 - alpha) / beta))
                if I_mtr[i, j] == 1:
                    rc2.add_clause([Y[i, j]], weight=math.log((1 - beta) / alpha))

    if kmax > 0:
        for combo in combinations(range(num_mutations), kmax + 1):
            tmp = []
            for com in combo:
                tmp.append(-K[com])
            rc2.add_clause(tmp)

        for j in range(num_mutations):
            j

    delta

    s_time = time.time()
    variables = rc2.compute()
    e_time = time.time()
    running_time = e_time - s_time

    sol_Y = np.empty((num_cells, num_mutations), dtype=np.int8)
    numVar = 0
    for i in range(num_cells):
        for j in range(num_mutations):
            sol_Y[i, j] = variables[numVar] > 0
            numVar += 1

    removedMutsIDs = []
    if kmax > 0:
        sol_K = np.empty(num_mutations, dtype=np.int8)
        numVar = numVarY + numVarB + numVarZ
        for j in range(num_mutations):
            sol_K[j] = variables[numVar] > 0
            numVar += 1
            if sol_K[j]:
                removedMutsIDs.append(snvs[j])

    df_output = pd.DataFrame(sol_Y)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    df_output[removedMutsIDs] = 0

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output


def phiscsi_bulk(
    df_input,
    alpha,
    beta,
    kmax=0,
    vaf_info=None,
    delta=0.2,
    time_limit=86400,
    n_threads=1,
):
    """Solving using PhISCS-I (in single-cell mode with bulk and mutation elimination).

    a combinatorial approach for subperfect tumor phylogeny reconstruction
    via integrative use of single-cell and bulk sequencing data :cite:`PhISCS`.

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
    kmax : :obj:`int`, optional
        Max number of mutations to be eliminated, by default 0
    vaf_info : :class:`pandas.DataFrame`, optional
        Information about the variant allele frequency in bulk data
        The size is n_SNVs x n_samples, by default None
    delta : :obj:`float`, optional
        Delta parameter accounting for VAF variance, by default 0.2
    time_limit : :obj:`int`, optional
        Time limit of the Gurobi running in seconds, by default 86400 (one day)
    n_threads : :obj:`int`, optional
        Number of threads for Gurobi solver, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    gp, gp_is_not_imported = tsc.ul.import_gurobi()
    if gp_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    tsc.logg.info(
        f"running PhISCS-I-orig with alpha={alpha}, beta={beta}, kmax={kmax}, "
        f"vaf_info={vaf_info}, delta={delta}, time_limit={time_limit}, "
        f"n_threads={n_threads}"
    )

    cells = list(df_input.index)
    snvs = list(df_input.columns)
    df_input = df_input.replace("?", 3)
    df_input = df_input.astype(int)
    I_mtr = df_input.values

    model = gp.Model("ILP")
    model.Params.OutputFlag = 0
    model.Params.LogFile = ""
    model.Params.Threads = n_threads
    model.Params.TimeLimit = time_limit

    numCells = len(cells)
    numMutations = len(snvs)
    if vaf_info is not None:
        vaf_info = vaf_info.loc[snvs]
        sampleIDs = vaf_info.columns
        vaf_info.loc["NULL"] = vaf_info.shape[1] * [1]

    # Matrix Y is matrix of corrected (i.e. true) genotypes w.r.t. input SC matrix I
    Y = {}
    for c in range(numCells):
        for m in range(numMutations):
            Y[c, m] = model.addVar(vtype=gp.GRB.BINARY, name=f"Y({c},{m})")

    # Variables B control the existence of conflict between columns
    B = {}
    for p in range(numMutations + 1):
        for q in range(numMutations + 1):
            B[p, q, 1, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,1]"
            )
            B[p, q, 1, 0] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,0]"
            )
            B[p, q, 0, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},0,1]"
            )

    K = {}
    for m in range(numMutations + 1):
        K[m] = model.addVar(vtype=gp.GRB.BINARY, name=f"K[{m}]")
    model.addConstr(K[numMutations] == 0)  # null mutation can not be eliminated

    # A[p,q] = 1 if p is ancestor of q
    A = {}
    if vaf_info is not None:
        for p in range(
            numMutations + 1
        ):  # mutation with index numMutation is null mutation
            for q in range(numMutations + 1):
                A[p, q] = model.addVar(vtype=gp.GRB.BINARY, obj=0, name=f"A[{p},{q}]")

    model.update()

    # number of eliminated columns is upper bounded by user provided constant
    model.addConstr(gp.quicksum(K[m] for m in range(numMutations)) <= kmax)

    # Enforce three gametes rule
    for i in range(numCells):
        for p in range(numMutations):
            for q in range(numMutations):
                model.addConstr(Y[i, p] + Y[i, q] - B[p, q, 1, 1] <= 1)
                model.addConstr(-Y[i, p] + Y[i, q] - B[p, q, 0, 1] <= 0)
                model.addConstr(Y[i, p] - Y[i, q] - B[p, q, 1, 0] <= 0)

    # Null mutation present in each cell
    for p in range(numMutations + 1):
        model.addConstr(B[p, numMutations, 1, 0] == 0)

    # Forbid conflict between columns (three gametes rule)
    for p in range(numMutations):
        for q in range(numMutations):
            model.addConstr(
                B[p, q, 0, 1] + B[p, q, 1, 0] + B[p, q, 1, 1] <= 2 + K[p] + K[q]
            )

    # Constraints for integrating VAF obtained from bulk data into the model
    if vaf_info is not None:
        for p in range(numMutations):
            for q in range(p + 1, numMutations):
                model.addConstr(A[p, q] + A[q, p] <= 1 - K[p])
                model.addConstr(A[p, q] + A[q, p] <= 1 - K[q])
                model.addConstr(A[p, q] + A[q, p] >= B[p, q, 1, 1] - K[p] - K[q])
        for p in range(numMutations + 1):
            model.addConstr(A[p, p] == 0)
            for q in range(numMutations + 1):
                model.addConstr(A[p, q] <= 1 - K[p])
                model.addConstr(A[p, q] <= 1 - K[q])

                if p < q:
                    model.addConstr(A[p, q] + A[q, p] <= 1)

                model.addConstr(A[p, q] + B[p, q, 0, 1] <= 1 + K[p] + K[q])
                model.addConstr(
                    B[p, q, 1, 0] + B[p, q, 1, 1] - A[p, q] <= 1 + K[p] + K[q]
                )

                for sampleID in sampleIDs:
                    VAF_p = float(vaf_info.iloc[p][sampleID])
                    VAF_q = float(vaf_info.iloc[q][sampleID])
                    model.addConstr(A[p, q] * VAF_p * (1 + delta) >= A[p, q] * VAF_q)

                    for r in range(numMutations + 1):
                        if r == q:
                            continue
                        VAF_r = float(vaf_info.iloc[r][sampleID])
                        # Constraint 2
                        model.addConstr(
                            VAF_p * (1 + delta)
                            >= VAF_q * (A[p, q] - A[r, q] - A[q, r])
                            + VAF_r * (A[p, r] - A[r, q] - A[q, r])
                        )

                for r in range(numMutations + 1):
                    if r == q:
                        continue
                    # Constraint 1.d
                    model.addConstr(A[p, r] >= A[p, q] + A[q, r] - 1)

            candidateAncestors = list(range(numMutations + 1))
            candidateAncestors.remove(p)

            if p < numMutations:
                model.addConstr(
                    gp.quicksum(A[s, p] for s in candidateAncestors) >= 1 - K[p]
                )
            elif p == numMutations:
                model.addConstr(gp.quicksum(A[s, p] for s in candidateAncestors) == 0)
            else:
                tsc.logg.error("p index is out of range")

    # Defining the objective function
    objective = 0
    for j in range(numMutations):
        numZeros = 0
        numOnes = 0
        for i in range(numCells):
            if I_mtr[i][j] == 0:
                numZeros += 1
                objective += np.log(beta / (1 - alpha)) * Y[i, j]
            elif I_mtr[i][j] == 1:
                numOnes += 1
                objective += np.log((1 - beta) / alpha) * Y[i, j]

        objective += numZeros * np.log(1 - alpha)
        objective += numOnes * np.log(alpha)
        objective -= K[j] * (
            numZeros * np.log(1 - alpha)
            + numOnes * (np.log(alpha) + np.log((1 - beta) / alpha))
        )

    model.setObjective(objective, gp.GRB.MAXIMIZE)

    s_time = time.time()
    model.optimize()
    e_time = time.time()
    running_time = e_time - s_time

    removedMutsIDs = []
    sol_K = []
    for j in range(numMutations):
        sol_K.append(K[j].X > 0.5)
        if sol_K[j] == 1:
            removedMutsIDs.append(snvs[j])

    sol_Y = np.zeros((numCells, numMutations), dtype=np.int8)
    for i in range(numCells):
        for j in range(numMutations):
            sol_Y[i, j] = Y[i, j].X > 0.5

    df_output = pd.DataFrame(sol_Y)
    df_output.columns = snvs
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    df_output[removedMutsIDs] = 0
    # df_output.drop(removedMutsIDs, axis=1, inplace=True)

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output


def phiscs_readcount(adata, alpha, beta, time_limit=86400, n_threads=1):
    """Solving using PhISCS-ReadCount (only in single-cell mode, no bulk).

    a combinatorial approach for subperfect tumor phylogeny reconstruction
    via integrative use of single-cell and bulk sequencing data :cite:`PhISCS`.

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
    time_limit : :obj:`int`, optional
        Time limit of the Gurobi running in seconds, by default 86400 (one day)
    n_threads : :obj:`int`, optional
        Number of threads for Gurobi solver, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """

    gp, gp_is_not_imported = tsc.ul.import_gurobi()
    if gp_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    tsc.logg.info(
        f"running PhISCS-readcout with alpha={alpha}, beta={beta},"
        f" time_limit={time_limit}"
    )

    PROB_SEQ_ERROR = 0.001
    BETABINOM_ALPHA = 1.0
    BETABINOM_BETA = 1.0

    def _prob_absent(v, t):
        prob = binom.pmf(v, t, PROB_SEQ_ERROR)
        return prob

    def _prob_present(v, t):
        prob = pmf_BetaBinomial(v, t, BETABINOM_ALPHA, BETABINOM_BETA)
        return prob

    df_input = adata.to_df()
    cells = adata.obs_names
    muts = adata.var_names
    T = adata.layers["total"]
    V = adata.layers["mutant"]

    B_mu = np.zeros((len(cells), len(muts)), dtype=np.float64)
    B_delta = np.zeros((len(cells), len(muts)), dtype=np.float64)
    for i in range(len(cells)):
        for j in range(len(muts)):
            if T[i, j] != 0:
                B_mu[i, j] = _prob_absent(V[i, j], T[i, j])
                B_delta[i, j] = _prob_present(V[i, j], T[i, j])

    model = gp.Model("ILP")
    model.Params.OutputFlag = 0
    model.Params.LogFile = ""
    model.Params.Threads = n_threads
    model.Params.TimeLimit = time_limit
    Y = {}
    B = {}
    for c in range(len(cells)):
        for m in range(len(muts)):
            Y[c, m] = model.addVar(vtype=gp.GRB.BINARY, name=f"Y({c},{m})")
    for p in range(len(muts)):
        for q in range(len(muts)):
            B[p, q, 1, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,1]"
            )
            B[p, q, 1, 0] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,0]"
            )
            B[p, q, 0, 1] = model.addVar(
                vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},0,1]"
            )
    for p in range(len(muts)):
        for q in range(len(muts)):
            model.addConstr(B[p, q, 0, 1] + B[p, q, 1, 0] + B[p, q, 1, 1] <= 2)
            for i in range(len(cells)):
                model.addConstr(Y[i, p] + Y[i, q] - B[p, q, 1, 1] <= 1)
                model.addConstr(-Y[i, p] + Y[i, q] - B[p, q, 0, 1] <= 0)
                model.addConstr(Y[i, p] - Y[i, q] - B[p, q, 1, 0] <= 0)
    objective = 0
    for j in range(len(muts)):
        for i in range(len(cells)):
            if T[i, j] != 0:
                objective += (1 - Y[i, j]) * np.log(
                    B_mu[i, j] * (1 - alpha) + B_delta[i, j] * (alpha)
                )
                objective += (Y[i, j]) * np.log(
                    B_mu[i, j] * (beta) + B_delta[i, j] * (1 - beta)
                )
    model.setObjective(objective, gp.GRB.MAXIMIZE)

    s_time = time.time()
    model.optimize()
    e_time = time.time()
    running_time = e_time - s_time

    sol_Y = np.zeros((len(cells), len(muts)), dtype=np.int8)
    for i in range(len(cells)):
        for j in range(len(muts)):
            sol_Y[i, j] = Y[i, j].X > 0.5

    df_output = pd.DataFrame(sol_Y)
    df_output.columns = muts
    df_output.index = cells
    df_output.index.name = "cellIDxmutID"

    tsc.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
