import itertools
import math

import trisicell as tsc
from trisicell.external.gpps._utils_ilp import expand_name

__author__ = "Simone Ciccolella"
__date__ = "11/30/21"


def gpps_ilp(
    input_matrix,
    alpha,
    beta,
    k_dollo,
    max_del=-1,
    time_limit=86400,
    n_threads=1,
):
    gp, gp_is_not_imported = tsc.ul.import_gurobi()
    if gp_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    # ==================================================================#
    # ========================= PREREQUISITES ==========================#
    # ==================================================================#

    # parser.add_argument(
    #     "--mps",
    #     action="store_true",
    #     help="This will output the model in MPS format instead of running the solver",
    # )

    # ----------------------Initialize program----------------------#
    # Fixed parameters
    # num_samples = len(input_matrix)
    # num_clones = int(num_mutations * n_clones)
    num_clones = len(input_matrix)
    num_mutations = len(input_matrix[0])

    # tsc.logg.debug("Num samples: %d" % num_samples)
    tsc.logg.debug("Num mutations: %d" % num_mutations)
    tsc.logg.debug("Num clones: %d" % num_clones)

    # ==================================================================#
    # ========================== GUROBI MODEL ==========================#
    # ==================================================================#
    model = gp.Model("Parsimony Phylogeny Model")
    model.Params.OutputFlag = 0
    model.Params.LogFile = ""
    model.setParam("Threads", n_threads)
    model.setParam("TimeLimit", time_limit)

    # ---------------------------------------------------#
    # ------------------- VARIABLES ---------------------#
    # ---------------------------------------------------#

    # -----------Variable Y and B---------------

    lalpha = math.log(alpha)
    lbeta = math.log(beta)
    l_alpha = math.log(1 - alpha)
    l_beta = math.log(1 - beta)

    tsc.logg.debug(lalpha, lbeta, l_alpha, l_beta)

    tsc.logg.debug("Generating variables I, F, w, P.")
    I_mtx = {}
    F = {}
    fminus = {}
    P = {}
    # False positive should be only a few
    FP = {}

    for row_index, row in enumerate(input_matrix):
        I_mtx[row_index] = {}
        F[row_index] = {}
        fminus[row_index] = {}
        P[row_index] = {}
        FP[row_index] = {}
        for col_index, cell in enumerate(row):
            # P
            names = expand_name(str(col_index), 1, k_dollo)
            for name in names:
                P[row_index][name] = model.addVar(
                    vtype=gp.GRB.BINARY, obj=0, name=f"P{row_index}-{name}"
                )
            if cell < 2:
                # I
                I_mtx[row_index][col_index] = model.addVar(
                    vtype=gp.GRB.BINARY,
                    obj=0,
                    name=f"I{row_index}-{col_index}",
                )
                model.update()
                model.addConstr(
                    I_mtx[row_index][col_index] == cell,
                    f"constr_I{row_index}-{col_index}",
                )
                # F and fminus. fminus is equal to 1-F
                if cell == 0:
                    F[row_index][col_index] = model.addVar(
                        vtype=gp.GRB.BINARY,
                        obj=lalpha,
                        name=f"F{row_index}-{col_index}",
                    )
                    fminus[row_index][col_index] = model.addVar(
                        vtype=gp.GRB.BINARY,
                        obj=l_beta,
                        name=f"f{row_index}-{col_index}",
                    )
                if cell == 1:
                    F[row_index][col_index] = model.addVar(
                        vtype=gp.GRB.BINARY,
                        obj=l_alpha,
                        name=f"F{row_index}-{col_index}",
                    )
                    fminus[row_index][col_index] = model.addVar(
                        vtype=gp.GRB.BINARY,
                        obj=lbeta,
                        name=f"f{row_index}-{col_index}",
                    )
                    FP[row_index][col_index] = model.addVar(
                        vtype=gp.GRB.BINARY,
                        obj=0,
                        name=f"FP{row_index}-{col_index}",
                    )
                    model.update()
                    model.addConstr(
                        FP[row_index][col_index] == 1 - P[row_index][names[0]],
                        name=f"constr_FP{row_index}-{col_index}",
                    )
                model.update()
                model.addConstr(
                    F[row_index][col_index] == 1 - fminus[row_index][col_index],
                    name=f"constr_def_fminus{row_index}-{col_index}",
                )
                model.addConstr(
                    F[row_index][col_index]
                    == P[row_index][names[0]]
                    - gp.quicksum(P[row_index][name] for name in names[1:]),
                    name=f"constr_balance_F{row_index}-{col_index}",
                )
                model.addConstr(
                    P[row_index][names[0]]
                    >= gp.quicksum(P[row_index][name] for name in names[1:]),
                    name=f"constr_imbalance_P{row_index}-{col_index}",
                )

    # There are only a few false positives
    model.addConstr(
        4 >= gp.quicksum(FP[r][c] for r in FP.keys() for c in FP[r].keys()),
        name="constr_few_FP",
    )

    model.update()

    tsc.logg.debug("Generating variables B.")
    B = {}
    columns = list(
        itertools.chain.from_iterable(
            [expand_name(str(name), 1, k_dollo) for name in range(num_mutations)]
        )
    )
    for p in columns:
        B[p] = {}
    for p, q in itertools.combinations(columns, 2):
        B[p][q] = {}
        B[p][q]["01"] = model.addVar(vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},0,1]")
        B[p][q]["10"] = model.addVar(vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,0]")
        B[p][q]["11"] = model.addVar(vtype=gp.GRB.BINARY, obj=0, name=f"B[{p},{q},1,1]")
        model.update()
        for row_index, _ in enumerate(input_matrix):
            model.addConstr(
                B[p][q]["01"] >= P[row_index][q] - P[row_index][p],
                f"constr_B01-{p}-{q}-{row_index}",
            )
            model.addConstr(
                B[p][q]["10"] >= P[row_index][p] - P[row_index][q],
                f"constr_B10-{p}-{q}-{row_index}",
            )
            model.addConstr(
                B[p][q]["11"] >= P[row_index][p] + P[row_index][q] - 1,
                f"constr_B11-{p}-{q}-{row_index}",
            )
        model.addConstr(
            B[p][q]["01"] + B[p][q]["10"] + B[p][q]["11"] <= 2,
            f"constr_sum_B{p}-{q}",
        )

    deletions = {}
    del_names = []
    for p in columns:
        if "-" in p:
            del_names.append(p)
            deletions[p] = model.addVar(vtype=gp.GRB.BINARY, obj=0, name=f"Del[{p}]")
            for row_index, _ in enumerate(input_matrix):
                model.addConstr(deletions[p] >= P[row_index][p])

    if max_del != -1:
        model.addConstr(gp.quicksum(deletions[x] for x in del_names) <= max_del)

    model.update()
    model.modelSense = gp.GRB.MAXIMIZE
    model.update()

    # ---------------------------------------------------#
    # -------------------- OPTIMIZE ---------------------#
    # ---------------------------------------------------#
    tsc.logg.debug("#----- GUROBI OPTIMIZATION ----#")
    model.optimize()

    # ==================================================================#
    # ======================= POST OPTIMIZATION ========================#
    # ==================================================================#

    output_matrix = []
    for row_index, row in enumerate(input_matrix):
        row_out = []
        for col_index, _ in enumerate(row):
            names = expand_name(str(col_index), 1, k_dollo)
            for name in names:
                row_out.append(int(float(P[row_index][name].X)))
        output_matrix.append(row_out)

    return output_matrix
