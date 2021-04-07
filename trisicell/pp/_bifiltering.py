import math

import numpy as np

import trisicell as tsc


def bifiltering(df, cellr, mutr, time_limit=3600):
    """Bi-filtering to find maximally inforemed submatrix.

    This function runs an ILP to find maximally inforemed submatrix
    where the number of mutant genotypes is maximized.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        The input noisy genotype matrix where entries are 0,1 and 3.
    cellr : :obj:`float`
        ratio for picking how many cells.
    mutr : :obj:`float`
        ratio for picking how many mutations.
    time_limit : :obj:`int`, optional
        Time limit for the ILP solver, by default 3600.

    Returns
    -------
    :class:`pandas.DataFrame`
        The output gentoype submatrix.
    """

    gp, gp_is_not_imported = tsc.ul.import_gurobi()
    if gp_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    M = df.values.copy()
    #### only mutant
    # M[M != 1] = -1

    #### mutant and reference
    M[M != 3] = 1
    M[M == 3] = -1
    n, m = M.shape
    n_cells, n_sites = math.floor(cellr * n), math.floor(mutr * m)

    model = gp.Model("bi-filtering")
    model.Params.OutputFlag = 0
    model.Params.LogFile = ""
    model.Params.Threads = 1
    model.Params.TimeLimit = time_limit

    tsc.logg.info("1. declaring variables", time=True)
    A, C, S = {}, {}, {}
    for i in range(n):
        A[i] = {}
        for j in range(m):
            A[i][j] = model.addVar(vtype=gp.GRB.BINARY)
    for i in range(n):
        C[i] = model.addVar(vtype=gp.GRB.BINARY)
        for j in range(m):
            S[j] = model.addVar(vtype=gp.GRB.BINARY)
    # A = model.addMVar((n, m), vtype=gp.GRB.BINARY)
    # C = model.addMVar(n, vtype=gp.GRB.BINARY)
    # S = model.addMVar(m, vtype=gp.GRB.BINARY)

    tsc.logg.info("2. adding constraints", time=True)
    for i in range(n):
        for j in range(m):
            model.addConstr(A[i][j] <= C[i])
            model.addConstr(A[i][j] <= S[j])
            model.addConstr(C[i] + S[j] - 1 <= A[i][j])
    # for i in range(n):
    #     for j in range(m):
    #         model.addConstr(A[i, j] <= C[i])
    #         model.addConstr(A[i, j] <= S[j])
    #         model.addConstr(C[i] + S[j] - 1 <= A[i, j])

    model.addConstr(gp.quicksum(C[i] for i in range(n)) == n_cells)
    model.addConstr(gp.quicksum(S[j] for j in range(m)) == n_sites)
    # model.addConstr(sum(C[i] for i in range(n)) == n_cells)
    # model.addConstr(sum(S[j] for j in range(m)) == n_sites)

    model.update()

    tsc.logg.info("3. adding objective", time=True)
    model.setObjective(
        gp.quicksum(A[i][j] * M[i, j] for i in range(n) for j in range(m)),
        gp.GRB.MAXIMIZE,
    )

    tsc.logg.info("4. start solving", time=True)
    model.optimize()

    a = np.array([i for i in range(n) if C[i].X > 0])
    b = np.array([j for j in range(m) if S[j].X > 0])
    cells = df.index[a]
    sites = df.columns[b]

    return df.loc[cells, sites].copy()
