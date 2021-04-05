import copy
import itertools
import time

import numpy as np
import pandas as pd
import pybnb
import scipy.sparse as sp
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF

import trisicell as tsc

rec_num = 0


def bnb(df_input, bounding):
    mpi4py, mpi4py_is_not_imporeted = tsc.ul.import_mpi4py()
    if mpi4py_is_not_imporeted:
        raise RuntimeError("Unable to import a package!")

    tsc.logg.info(f"running BnB with bounding={bounding}")
    matrix_input = df_input.values
    matrix_output = matrix_input.copy()
    na_value = 3
    bounding_algs = {
        "real": TwoSatBounding(
            heuristic_setting=None,
            n_levels=2,
            compact_formulation=False,
            na_value=na_value,
        ),  # Real Data
        "simulated": TwoSatBounding(
            heuristic_setting=[True, True, False, True, True],
            n_levels=1,
            compact_formulation=True,
            na_value=na_value,
        ),  # Simulation
    }

    bounding_alg = bounding_algs[bounding]

    s_time = time.time()
    flips = solve_by_BnB(matrix_input, na_value, bounding_alg)
    e_time = time.time()
    running_time = e_time - s_time

    for k in flips:
        matrix_output[k] = 1
    matrix_output[np.where(matrix_output == na_value)] = 0

    df_output = pd.DataFrame(matrix_output)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    tsc.ul.stat(df_input, df_output, 0.1, 0.000000001, running_time)

    return df_output


def solve_by_BnB(matrix_in, na_value, bounding_alg):
    result = bnb_solve(matrix_in, bounding_algorithm=bounding_alg, na_value=na_value)
    matrix_output = result[0]
    flips = []
    zero_one_flips = np.where((matrix_in != matrix_output) & (matrix_in != na_value))
    for i in range(len(zero_one_flips[0])):
        flips.append((zero_one_flips[0][i], zero_one_flips[1][i]))
    na_one_flips = np.where((matrix_output == 1) & (matrix_in == na_value))
    for i in range(len(na_one_flips[0])):
        flips.append((na_one_flips[0][i], na_one_flips[1][i]))
    return flips


def all_None(*args):
    return args.count(None) == len(args)


def calculate_column_intersections(matrix, for_loop=False, row_by_row=False):
    ret = np.empty((matrix.shape[1], matrix.shape[1]), dtype=np.bool)
    mask_1 = matrix == 1

    if for_loop:
        for p in range(matrix.shape[1]):
            # even though the diagonals are not necessary, I keep it for ease of debugging
            for q in range(p, matrix.shape[1]):
                ret[p, q] = np.any(np.logical_and(mask_1[:, p], mask_1[:, q]))
                ret[q, p] = ret[p, q]
    elif row_by_row:
        ret[:, :] = 0
        for r in range(matrix.shape[0]):
            one_columns = mask_1[r]
            ret[np.ix_(one_columns, one_columns)] = True
    return ret


def zero_or_na(vec, na_value=-1):
    return np.logical_or(vec == 0, vec == na_value)


def make_sure_variable_exists(
    memory_matrix, row, col, num_var_F, map_f2ij, var_list, na_value
):
    if memory_matrix[row, col] < 0:
        num_var_F += 1
        map_f2ij[num_var_F] = (row, col)
        memory_matrix[row, col] = num_var_F
        var_list.append(num_var_F)
    return num_var_F


def get_effective_matrix(I, delta01, delta_na_to_1, change_na_to_0=False):
    x = np.array(I + delta01, dtype=np.int8)
    if delta_na_to_1 is not None:
        na_indices = delta_na_to_1.nonzero()
        x[
            na_indices
        ] = 1  # should have been (but does not accept): x[na_indices] = delta_na_to_1[na_indices]
    if change_na_to_0:
        x[np.logical_and(x != 0, x != 1)] = 0
    return x


def make_twosat_model_from_np(
    constraints,
    F,
    zero_vars,
    na_vars,
    eps=None,
    heuristic_setting=None,
    compact_formulation=True,
):

    if eps is None:
        eps = 1 / (len(zero_vars) + len(na_vars))

    if heuristic_setting is None:
        rc2 = RC2(WCNF())
    else:
        assert len(heuristic_setting) == 5
        rc2 = RC2(
            WCNF(),
            adapt=heuristic_setting[0],
            exhaust=heuristic_setting[1],
            incr=heuristic_setting[2],
            minz=heuristic_setting[3],
            trim=heuristic_setting[4],
        )

    if not compact_formulation:
        # hard constraints Z_a,p or Z_b,q
        for constr_ind in range(constraints[0].shape[0]):
            constraint = constraints[0][constr_ind]
            a, p, b, q = constraint.flat
            # print(constraint, F.shape)
            # print(a, p, b, q)
            rc2.add_clause([F[a, p], F[b, q]])
        if len(constraints) >= 2:
            # hard constraints Z_a,p or Z_b,q or -Z_c,d
            for constr_ind in range(constraints[1].shape[0]):
                constraint = constraints[1][constr_ind]
                a, p, b, q, c, d = constraint.flat
                # print(a, p, b, q, c, d)
                rc2.add_clause([F[a, p], F[b, q], -F[c, d]])
    else:
        # hard constraints Z_a,p or (sign) b_pq
        for constr_ind in range(constraints[0].shape[0]):
            constraint = constraints[0][constr_ind]
            row, col, b_pq, sign = constraint.flat
            rc2.add_clause([F[row, col], sign * b_pq])
        if len(constraints) >= 2:
            # hard constraints Z_a,p or Z_b,q or -Z_c,d
            for constr_ind in range(constraints[1].shape[0]):
                constraint = constraints[1][constr_ind]
                row, col, c_pq0, c_pq1 = constraint.flat
                # if Z_rc is True at least one of p, q should become active
                # E.g., c_pq0 be False
                rc2.add_clause([-F[row, col], -c_pq0, -c_pq1])
                # if c_pq0 is False then Z_rc has to be flipped
                rc2.add_clause([F[row, col], c_pq0])

    # soft constraints for zero variables
    for var in zero_vars:
        rc2.add_clause([-var], weight=1)

    if eps > 0:
        # soft constraints for zero variables
        for var in na_vars:
            rc2.add_clause([-var], weight=eps)

    return rc2


def twosat_solver(
    matrix,
    cluster_rows=False,
    cluster_cols=False,
    only_descendant_rows=False,
    na_value=None,
    leave_nas_if_zero=False,
    return_lb=False,
    heuristic_setting=None,
    n_levels=2,
    eps=0,
    compact_formulation=False,
):
    global rec_num
    rec_num += 1
    assert not cluster_rows, "Not implemented yet"
    assert not cluster_cols, "Not implemented yet"
    assert not only_descendant_rows, "Not implemented yet"
    model_time = 0
    opt_time = 0
    start_time = time.time()

    return_value = make_constraints_np_matrix(
        matrix,
        n_levels=n_levels,
        na_value=na_value,
        compact_formulation=compact_formulation,
    )
    model_time += time.time() - start_time
    F, map_f2ij, zero_vars, na_vars, hard_constraints, col_pair = (
        return_value.F,
        return_value.map_f2ij,
        return_value.zero_vars,
        return_value.na_vars,
        return_value.hard_constraints,
        return_value.col_pair,
    )

    if col_pair is not None:
        icf = False
    elif return_value.complete_version:
        icf = True
    else:
        icf = None

    final_output = None
    lower_bound = 0
    if icf:
        final_output, total_time = matrix.copy(), 0
    else:
        start_time = time.time()
        rc2 = make_twosat_model_from_np(
            hard_constraints,
            F,
            zero_vars,
            na_vars,
            eps,
            heuristic_setting,
            compact_formulation=compact_formulation,
        )
        model_time += time.time() - start_time

        a = time.time()
        variables = rc2.compute()
        b = time.time()
        opt_time += b - a
        output_matrix = matrix.copy()
        output_matrix = output_matrix.astype(np.int8)

        for var_ind in range(len(variables)):
            if (
                0 < variables[var_ind] and variables[var_ind] in map_f2ij
            ):  # if 0 or 2 make it one
                output_matrix[map_f2ij[variables[var_ind]]] = 1
                if matrix[map_f2ij[variables[var_ind]]] != na_value:
                    lower_bound += 1
        # I don't change 2s to 0s here keep them 2 for next time

        # For recursion I set off all sparsification parameters
        # Also I want na->0 to stay na for the recursion regardless of original input for leave_nas_if_zero
        # I am also not passing eps here to wrap up the recursion soon

        Orec, rec_model_time, rec_opt_time = twosat_solver(
            output_matrix,
            na_value=na_value,
            heuristic_setting=None,
            n_levels=n_levels,
            leave_nas_if_zero=True,
            compact_formulation=compact_formulation,
        )
        model_time += rec_model_time
        opt_time += rec_opt_time

        if not leave_nas_if_zero:
            Orec[Orec == na_value] = 0
        final_output = Orec

    if return_lb:
        return final_output, model_time, opt_time, lower_bound
    else:
        return final_output, model_time, opt_time


def make_constraints_np_matrix(
    matrix,
    constraints=None,
    n_levels=2,
    na_value=None,
    row_coloring=None,
    col_coloring=None,
    probability_threshold=None,
    fn_rate=None,
    column_intersection=None,
    compact_formulation=True,
):
    """
    Returns a "C x 2 x 2" matrix where C is the number of extracted constraints each constraints is of the form:
    ((r1, c1), (r2, c2)) and correspond to Z_{r1, c1} or Z{r2, c2}
    :param matrix: A binary matrix cellsXmutations
    :param constraints: If not None instead of evaluating the whole matrix it will only look at potential constraints
    :param level: The type of constraints to add
    :param na_value:
    :param row_coloring: Only constraints that has the same row coloring will be used
    :param col_coloring: Only constraints that has the same column coloring will be used
    :param probability_threshold:
    :param fn_rate:
    :return:
    """
    # TD: Take decendence analysis out of here?
    # TD: how to reuse constraints input
    from collections import namedtuple

    assert (probability_threshold is None) == (fn_rate is None)
    descendance_analysis = probability_threshold is not None
    assert 1 <= n_levels <= 2, "not implemented yet"

    # means none of scarification ideas have been used
    complete_version = all_None(
        row_coloring, col_coloring, probability_threshold, fn_rate
    )

    soft_cnst_num = 0
    hard_constraints = [[] for _ in range(n_levels)]  # an empty list each level
    if descendance_analysis:
        # dictionary for lazy calculation of decadence:
        descendent_dict = dict()

    # variables for each zero
    F = -np.ones(matrix.shape, dtype=np.int64)
    num_var_F = 0
    map_f2ij = dict()
    zero_vars = list()
    na_vars = list()
    if compact_formulation:
        B_vars_offset = matrix.shape[0] * matrix.shape[1] + 1
        num_var_B = 0
        map_b2ij = dict()
        if n_levels >= 2:
            C_vars_offset = B_vars_offset + matrix.shape[1] * matrix.shape[1] + 1
            num_var_C = 0
            map_c2ij = dict()

    col_pair = None
    pair_cost = 0

    if column_intersection is None:
        column_intersection = calculate_column_intersections(matrix, row_by_row=True)
        # column_intersection = calculate_column_intersections(matrix, for_loop=True)
    for p in range(matrix.shape[1]):
        for q in range(p + 1, matrix.shape[1]):
            if column_intersection[p, q]:  # p and q has intersection
                # TD: check col_coloring here
                r01 = np.nonzero(
                    np.logical_and(
                        zero_or_na(matrix[:, p], na_value=na_value), matrix[:, q] == 1
                    )
                )[0]
                r10 = np.nonzero(
                    np.logical_and(
                        matrix[:, p] == 1, zero_or_na(matrix[:, q], na_value=na_value)
                    )
                )[0]
                cost = min(len(r01), len(r10))
                if cost > pair_cost:  # keep best pair to return as auxiliary info
                    # print("------------", cost, (p, q), len(r01), len(r10), column_intersection[p, q])
                    col_pair = (p, q)
                    pair_cost = cost
                if cost > 0:  # don't do anything if one of r01 or r10 is empty
                    if (
                        not compact_formulation
                    ):  # len(r01) * len(r10) many constraints will be added
                        for a, b in itertools.product(r01, r10):
                            # TD: check row_coloring
                            for row, col in [
                                (a, p),
                                (b, q),
                            ]:  # make sure the variables for this are made
                                var_list = (
                                    zero_vars if matrix[row, col] == 0 else na_vars
                                )
                                num_var_F = make_sure_variable_exists(
                                    F, row, col, num_var_F, map_f2ij, var_list, na_value
                                )
                            hard_constraints[0].append(
                                [[a, p], [b, q]]
                            )  # at least one of them should be flipped
                    else:  # compact formulation: (r01 + r10) number of new constraints will be added
                        # define new B variable
                        b_pq = B_vars_offset + num_var_B
                        num_var_B += 1
                        for row_list, col, sign in zip((r01, r10), (p, q), (1, -1)):
                            for row in row_list:
                                var_list = (
                                    zero_vars if matrix[row, col] == 0 else na_vars
                                )
                                num_var_F = make_sure_variable_exists(
                                    F, row, col, num_var_F, map_f2ij, var_list, na_value
                                )
                                hard_constraints[0].append([row, col, b_pq, sign])
                                # this will be translated to (Z_ap or (sign)B_pq)
            elif n_levels >= 2:
                r01 = np.nonzero(
                    np.logical_and(
                        zero_or_na(matrix[:, p], na_value=na_value), matrix[:, q] == 1
                    )
                )[0]
                r10 = np.nonzero(
                    np.logical_and(
                        matrix[:, p] == 1, zero_or_na(matrix[:, q], na_value=na_value)
                    )
                )[0]
                cost = min(len(r01), len(r10))
                if cost > 0:  # don't do anything if one of r01 or r10 is empty
                    if not compact_formulation:
                        # len(r01) * len(r10) * (len(r01) * len(r10)) many constraints will be added
                        x = np.empty((r01.shape[0] + r10.shape[0], 2), dtype=np.int)
                        x[: len(r01), 0] = r01
                        x[: len(r01), 1] = p
                        x[-len(r10) :, 0] = r10
                        x[-len(r10) :, 1] = q

                        for a, b, ind in itertools.product(r01, r10, range(x.shape[0])):
                            for row, col in [
                                (a, p),
                                (b, q),
                                (x[ind, 0], x[ind, 1]),
                            ]:  # make sure the variables for this are made
                                # print(row, col)
                                var_list = (
                                    zero_vars if matrix[row, col] == 0 else na_vars
                                )
                                num_var_F = make_sure_variable_exists(
                                    F, row, col, num_var_F, map_f2ij, var_list, na_value
                                )
                            row = [[a, p], [b, q], [x[ind, 0], x[ind, 1]]]
                            if not np.array_equal(
                                row[0], row[2]
                            ) and not np.array_equal(row[1], row[2]):
                                hard_constraints[1].append(
                                    [[a, p], [b, q], [x[ind, 0], x[ind, 1]]]
                                )
                    else:  #  if compact_formulation: 2(r01 + r10) will be added
                        # define two new C variable
                        c_pq0 = C_vars_offset + num_var_C
                        num_var_C += 1
                        c_pq1 = C_vars_offset + num_var_C
                        num_var_C += 1
                        for row_list, col, sign in zip((r01, r10), (p, q), (1, -1)):
                            for row in row_list:
                                var_list = (
                                    zero_vars if matrix[row, col] == 0 else na_vars
                                )
                                num_var_F = make_sure_variable_exists(
                                    F, row, col, num_var_F, map_f2ij, var_list, na_value
                                )
                                if sign == 1:
                                    hard_constraints[1].append([row, col, c_pq0, c_pq1])
                                    # this will be translated to (~Z_ap or ~c_pq0 or ~c_pq1)
                                    # and (Z_ap or c_pq0)
                                else:
                                    hard_constraints[1].append([row, col, c_pq1, c_pq0])
                                    # this will be translated to (~Z_ap or ~c_pq0 or ~c_pq1) (the same)
                                    # and (Z_ap or c_pq1) (different)

    # TD: when using this make sure to put an if to say if the model is small and
    return_type = namedtuple(
        "ReturnType",
        "F map_f2ij zero_vars na_vars hard_constraints col_pair complete_version",
    )
    for ind in range(n_levels):
        hard_constraints[ind] = np.array(hard_constraints[ind], dtype=np.int)
    return return_type(
        F, map_f2ij, zero_vars, na_vars, hard_constraints, col_pair, complete_version
    )


def is_conflict_free_gusfield_and_get_two_columns_in_coflicts(I, na_value):
    def sort_bin(a):
        b = np.transpose(a)
        b_view = np.ascontiguousarray(b).view(
            np.dtype((np.void, b.dtype.itemsize * b.shape[1]))
        )
        idx = np.argsort(b_view.ravel())[::-1]
        c = b[idx]
        return np.transpose(c), idx

    Ip = I.copy()
    Ip[Ip == na_value] = 0
    O, idx = sort_bin(Ip)
    # TD: delete duplicate columns
    # print(O, '\n')
    Lij = np.zeros(O.shape, dtype=int)
    for i in range(O.shape[0]):
        maxK = 0
        for j in range(O.shape[1]):
            if O[i, j] == 1:
                Lij[i, j] = maxK
                maxK = j + 1
    # print(Lij, '\n')
    Lj = np.amax(Lij, axis=0)
    # print(Lj, '\n')
    for i in range(O.shape[0]):
        for j in range(O.shape[1]):
            if O[i, j] == 1:
                if Lij[i, j] != Lj[j]:
                    return False, (idx[j], idx[Lj[j] - 1])
    return True, (None, None)


class BoundingAlgAbstract:
    def __init__(self):
        self.matrix = None
        self._extra_info = None
        self._extraInfo = {}
        self._times = {}
        self.na_support = False
        pass

    def reset(self, matrix):
        raise NotImplementedError("The method not implemented")

    def get_bound(self, delta):
        """
        This bound should include the flips done so far too
        delta: a sparse matrix with fliped ones
        """
        raise NotImplementedError("The method not implemented")

    def get_name(self):
        return type(self).__name__

    def get_state(self):
        return None

    def set_state(self, state):
        assert state is None
        pass

    def get_extra_info(self):
        """
        Some bounding algorithms can provide extra information after calling bounding.
        E.g.,
        return {"icf":True, "bestPair":(a,b)}
        """
        return copy.copy(self._extraInfo)

    def get_priority(self, till_here, this_step, after_here, icf=False):
        return -after_here

    def get_times(self):
        return self._times

    def get_init_node(self):
        return None


class TwoSatBounding(BoundingAlgAbstract):
    def __init__(
        self,
        priority_version=-1,
        cluster_rows=False,
        cluster_cols=False,
        only_descendant_rows=False,
        na_value=None,
        heuristic_setting=None,
        n_levels=2,
        eps=0,
        compact_formulation=False,
    ):
        """
        :param priority_version:
        """
        assert not cluster_rows, "Not implemented yet"
        assert not cluster_cols, "Not implemented yet"
        assert not only_descendant_rows, "Not implemented yet"

        self.priority_version = priority_version

        self.na_support = True
        self.na_value = na_value
        self.matrix = None
        self._times = None
        self.next_lb = None
        self.heuristic_setting = heuristic_setting
        self.n_levels = n_levels
        self.eps = eps  # only for upperbound
        self.compact_formulation = compact_formulation
        self.cluster_rows = cluster_rows
        self.cluster_cols = cluster_cols
        self.only_descendant_rows = only_descendant_rows

    def get_name(self):
        params = [
            type(self).__name__,
            self.priority_version,
            self.heuristic_setting,
            self.n_levels,
            self.eps,
            self.compact_formulation,
        ]
        params_str = map(str, params)
        return "_".join(params_str)

    def reset(self, matrix):
        self.matrix = matrix  # TD: make the model here and do small alterations later

        # self.na_value = infer_na_value(matrix)
        self._times = {"model_preparation_time": 0, "optimization_time": 0}

    def get_init_node(self):

        # def twosat_solver(matrix, cluster_rows=False, cluster_cols=False, only_descendant_rows=False,
        #                   na_value=None, leave_nas_if_zero=False, return_lb=False, heuristic_setting=None,
        #                   n_levels=2, eps=0, compact_formulation=True):
        #     pass

        node = pybnb.Node()
        solution, model_time, opt_time, lb = twosat_solver(
            self.matrix,
            cluster_rows=self.cluster_rows,
            cluster_cols=self.cluster_cols,
            only_descendant_rows=self.only_descendant_rows,
            na_value=self.na_value,
            leave_nas_if_zero=True,
            return_lb=True,
            heuristic_setting=None,
            n_levels=self.n_levels,
            eps=self.eps,
            compact_formulation=self.compact_formulation,
        )
        self._times["model_preparation_time"] += model_time
        self._times["optimization_time"] += opt_time

        nodedelta = sp.lil_matrix(np.logical_and(solution == 1, self.matrix == 0))
        node_na_delta = sp.lil_matrix(
            np.logical_and(solution == 1, self.matrix == self.na_value)
        )
        node.state = (
            nodedelta,
            True,
            None,
            nodedelta.count_nonzero(),
            self.get_state(),
            node_na_delta,
        )
        node.queue_priority = self.get_priority(
            till_here=-1, this_step=-1, after_here=-1, icf=True
        )
        self.next_lb = lb
        return node

    def get_bound(self, delta, delta_na=None):
        # make this dynamic when more nodes were getting explored
        if self.next_lb is not None:
            lb = self.next_lb
            self.next_lb = None
            return lb
        self._extraInfo = None
        current_matrix = get_effective_matrix(self.matrix, delta, delta_na)
        has_na = np.any(current_matrix == self.na_value)

        model_time = time.time()
        return_value = make_constraints_np_matrix(
            current_matrix,
            n_levels=self.n_levels,
            na_value=self.na_value,
            compact_formulation=self.compact_formulation,
        )
        F, map_f2ij, zero_vars, na_vars, hard_constraints, col_pair = (
            return_value.F,
            return_value.map_f2ij,
            return_value.zero_vars,
            return_value.na_vars,
            return_value.hard_constraints,
            return_value.col_pair,
        )

        if col_pair is not None:
            icf = False
        elif return_value.complete_version:
            icf = True
        else:
            icf = None  # not sure
        rc2 = make_twosat_model_from_np(
            hard_constraints,
            F,
            zero_vars,
            na_vars,
            eps=0,
            heuristic_setting=self.heuristic_setting,
            compact_formulation=self.compact_formulation,
        )

        model_time = time.time() - model_time
        self._times["model_preparation_time"] += model_time

        opt_time = time.time()
        variables = rc2.compute()
        opt_time = time.time() - opt_time
        self._times["optimization_time"] += opt_time

        result = 0
        for var_ind in range(len(variables)):
            if (
                variables[var_ind] > 0
                and abs(variables[var_ind]) in map_f2ij
                and self.matrix[map_f2ij[abs(variables[var_ind])]] == 0
            ):
                result += 1

        assert has_na or ((result == 0) == (col_pair is None)), f"{result}_{col_pair}"
        self._extraInfo = {
            "icf": icf,
            "one_pair_of_columns": col_pair,
        }
        ret = result + delta.count_nonzero()
        return ret

    def get_priority(self, till_here, this_step, after_here, icf=False):
        if icf:
            return self.matrix.shape[0] * self.matrix.shape[1] + 10
        else:
            sgn = np.sign(self.priority_version)
            pv_abs = self.priority_version * sgn
            if pv_abs == 1:
                return sgn * (till_here + this_step + after_here)
            elif pv_abs == 2:
                return sgn * (this_step + after_here)
            elif pv_abs == 3:
                return sgn * (after_here)
            elif pv_abs == 4:
                return sgn * (till_here + after_here)
            elif pv_abs == 5:
                return sgn * (till_here)
            elif pv_abs == 6:
                return sgn * (till_here + this_step)
            elif pv_abs == 7:
                return 0
        assert False, "get_priority did not return anything!"


class BnB(pybnb.Problem):
    def __init__(self, I, boundingAlg: BoundingAlgAbstract, na_value=None):
        self.na_value = na_value
        self.has_na = np.any(I == self.na_value)
        self.I = I
        self.delta = sp.lil_matrix(I.shape, dtype=np.int8)  # this can be coo_matrix too
        self.boundingAlg = boundingAlg
        self.delta_na = None
        if self.has_na:
            assert (
                boundingAlg.na_support
            ), "Input has N/A coordinates but bounding algorithm doesn't support it."
            self.delta_na = sp.lil_matrix(
                I.shape, dtype=np.int8
            )  # the coordinates with na that are decided to be 1
        (
            self.icf,
            self.colPair,
        ) = is_conflict_free_gusfield_and_get_two_columns_in_coflicts(self.I, na_value)
        self.boundingAlg.reset(I)
        self.node_to_add = self.boundingAlg.get_init_node()
        self.bound_value = self.boundingAlg.get_bound(self.delta)

    def sense(self):
        return pybnb.minimize

    def objective(self):
        if self.icf:
            return self.delta.count_nonzero()
        else:
            return pybnb.Problem.infeasible_objective(self)

    def bound(self):
        return self.bound_value

    def save_state(self, node):
        node.state = (
            self.delta,
            self.icf,
            self.colPair,
            self.bound_value,
            self.boundingAlg.get_state(),
            self.delta_na,
        )

    def load_state(self, node):
        (
            self.delta,
            self.icf,
            self.colPair,
            self.bound_value,
            boundingAlgState,
            self.delta_na,
        ) = node.state
        self.boundingAlg.set_state(boundingAlgState)

    def get_current_matrix(self):
        return get_effective_matrix(self.I, self.delta, self.delta_na)

    def branch(self):
        if self.icf:
            return

        need_for_new_nodes = True
        if self.node_to_add is not None:
            newnode = self.node_to_add
            self.node_to_add = None
            if (
                newnode.state[0].count_nonzero() == self.bound_value
            ):  # current_obj == lb => no need to explore
                need_for_new_nodes = False
            assert (
                newnode.queue_priority is not None
            ), "Right before adding a node its priority in the queue is not set!"
            yield newnode

        if need_for_new_nodes:
            p, q = self.colPair
            nf01 = None
            current_matrix = self.get_current_matrix()
            for col, colp in [(q, p), (p, q)]:
                node = pybnb.Node()
                nodedelta = copy.deepcopy(self.delta)
                node_na_delta = copy.deepcopy(self.delta_na)
                col1 = np.array(current_matrix[:, col], dtype=np.int8).reshape(-1)
                col2 = np.array(current_matrix[:, colp], dtype=np.int8).reshape(-1)
                rows01 = np.nonzero(np.logical_and(col1 == 0, col2 == 1))[0]
                rows21 = np.nonzero(np.logical_and(col1 == self.na_value, col2 == 1))[0]
                if (
                    len(rows01) + len(rows21) == 0
                ):  # nothing has changed! Dont add new node
                    continue
                nodedelta[rows01, col] = 1
                nf01 = nodedelta.count_nonzero()
                if self.has_na:
                    node_na_delta[rows21, col] = 1
                    new_bound = self.boundingAlg.get_bound(nodedelta, node_na_delta)
                else:
                    new_bound = self.boundingAlg.get_bound(nodedelta)

                node_icf, nodecol_pair = None, None
                extra_info = self.boundingAlg.get_extra_info()

                if extra_info is not None:
                    if "icf" in extra_info:
                        node_icf = extra_info["icf"]
                    if "one_pair_of_columns" in extra_info:
                        nodecol_pair = extra_info["one_pair_of_columns"]
                if node_icf is None:
                    x = get_effective_matrix(self.I, nodedelta, node_na_delta)
                    (
                        node_icf,
                        nodecol_pair,
                    ) = is_conflict_free_gusfield_and_get_two_columns_in_coflicts(
                        x, self.na_value
                    )

                node_bound_value = max(self.bound_value, new_bound)
                node.state = (
                    nodedelta,
                    node_icf,
                    nodecol_pair,
                    node_bound_value,
                    self.boundingAlg.get_state(),
                    node_na_delta,
                )
                node.queue_priority = self.boundingAlg.get_priority(
                    till_here=nf01 - len(rows01),
                    this_step=len(rows01),
                    after_here=new_bound - nf01,
                    icf=node_icf,
                )
                assert (
                    node.queue_priority is not None
                ), "Right before adding a node its priority in the queue is not set!"
                yield node


def bnb_solve(matrix, bounding_algorithm, na_value=None):
    problem1 = BnB(matrix, bounding_algorithm, na_value=na_value)
    solver = pybnb.solver.Solver()
    results1 = solver.solve(problem1, queue_strategy="custom", log=None)
    if results1.solution_status != "unknown":
        returned_delta = results1.best_node.state[0]
        returned_delta_na = results1.best_node.state[-1]
        returned_matrix = get_effective_matrix(
            matrix, returned_delta, returned_delta_na, change_na_to_0=True
        )
    else:
        returned_matrix = np.zeros((1, 1))
    # print("results1.nodes:  ", results1.nodes)
    return returned_matrix, results1.termination_condition
