import itertools
import time
from multiprocessing import Process, Queue

import numpy as np

__author__ = "Can Kizilkale"
__date__ = "3/19/21"

prfp = 5


def Reconstruct(
    df_input,
    Algchoice="FPNA",
    auto_tune=1,
    overlapp_coeff=0.3,
    hist_coeff=20,
    fnfp=100,
    fnc=1,
    postprcoef=5,
    post_fp=0.0001,
    post_fn=0.1,
    n_proc=7,
):
    global prfp
    prfp = postprcoef
    q = Queue()
    matrix_input = ReadFileNA(df_input.copy())
    matrix_input_raw = ReadFasis(df_input.copy())
    matrix_NA_est = Estimated_Matrix(df_input.copy())
    h_range = list(range(1, 100, 2))
    oc_range = [x / 10 for x in range(1, 5)]
    tune_var = itertools.product(h_range, oc_range)

    running_time = 0

    if Algchoice == "FN":
        s_time = time.time()
        matrix_recons = greedyPtreeNew(matrix_input.astype(bool))[1]
        e_time = time.time()
        running_time = e_time - s_time

    if Algchoice == "FPNA" and auto_tune == 0:
        s_time = time.time()
        apprx_ordr = sum(matrix_NA_est)
        matrix_recons = greedyPtreeNA(
            matrix_input.astype(bool), apprx_ordr, overlapp_coeff, hist_coeff
        )[2]
        e_time = time.time()
        running_time = e_time - s_time

    if Algchoice == "FPNA" and auto_tune == 1:
        tune_var = list(tune_var)
        proc_size = np.ceil(len(tune_var) / n_proc).astype(int)
        cpu_range = []

        for i in range(n_proc):
            s_i = i * proc_size
            e_i = s_i + proc_size
            if e_i > len(tune_var):
                e_i = len(tune_var)
            cpu_range.append(tune_var[s_i:e_i])
        s_time = time.time()
        p = [
            Process(
                target=Auto_fnfp,
                args=(
                    q,
                    cpu_range[i],
                    matrix_input,
                    matrix_NA_est,
                    matrix_input_raw,
                    fnfp,
                    fnc,
                ),
            )
            for _ in range(len(cpu_range))
        ]
        for i in p:
            i.start()
        for i in p:
            i.join()
        ret = []
        while not q.empty():
            ret.append(q.get())

        [m_r, d_min] = ret[0]
        matrix_recons = m_r
        for i in range(1, len(ret)):
            [m, d] = ret[i]
            if d < d_min:
                matrix_recons = m
                d_min = d

        e_time = time.time()
        running_time = e_time - s_time

    matrix_recons = postprocess_col(df_input, matrix_recons, pfn=post_fn, pfp=post_fp)

    return matrix_recons, running_time


def Auto_fnfp(q, tune_ran, m_input, m_NA_est, m_input_raw, fnfp, fnc):
    apprx_ordr = sum(m_NA_est)

    matrix_recon = greedyPtreeNA(
        m_input.astype(bool), apprx_ordr, tune_ran[0][1], tune_ran[0][0]
    )[2]
    matrix_rec_Temp = deleteNas(m_input_raw, matrix_recon)
    n_10 = np.sum(matrix_rec_Temp < m_input)
    n_01 = np.sum(matrix_rec_Temp > m_input)

    distance = fnfp * n_10 + fnc * n_01

    for [h_i, overlapp_coeff] in tune_ran:
        matrix_rec_i = greedyPtreeNA(
            m_input.astype(bool), apprx_ordr, overlapp_coeff, h_i
        )[2]
        matrix_rec_Temp = deleteNas(m_input_raw, matrix_rec_i)
        n_10 = np.sum(matrix_rec_Temp < m_input, dtype="int64")
        n_01 = np.sum(matrix_rec_Temp > m_input, dtype="int64")
        distance_i = fnfp * n_10 + fnc * n_01
        if distance_i < distance:
            matrix_recon = matrix_rec_i.copy()
            distance = distance_i
    q.put([matrix_recon, distance])


def deleteNas(M_in, M_out):
    M_o = M_out.copy()
    NA_position = np.argwhere(M_in == 3)
    for j in NA_position:
        M_o[j[0], j[1]] = 0
    return M_o


def greedyPtreeNew(M_input):
    M_copy = M_input.copy()

    ISet1 = np.argsort(sum(M_copy))
    ISet = []
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret = []
    while len(ISet) > 1:
        pivot_index = ISet[-1]
        Sremaining = ISet.copy()
        pivot_vector = M_copy[:, pivot_index]
        cum_vector = np.copy(pivot_vector)
        while_cont = 1

        while while_cont == 1:
            while_cont = 0

            for j in Sremaining:
                cap_j = cum_vector * M_copy[:, j]

                if np.any(cap_j):
                    cum_vector = cum_vector + M_copy[:, j]
                    while_cont = 1
                    Sremaining.remove(j)
        M_copy[:, pivot_index] = cum_vector
        ISet.remove(pivot_index)
        bret = np.argwhere(M_copy.astype(int) > M_input.astype(int))
    return [bret, M_copy]


def greedyPtreeNA(M_input, approx_order, oc, hc):
    overlap_coeff = oc
    hist_coeff = hc

    M_copy = M_input.copy()
    ISet1 = np.argsort(approx_order)
    ISet = []
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret = []
    pret = []

    while len(ISet) > 1:
        pivot_index = ISet[-1]
        Sremaining = ISet.copy()
        pivot_vector = M_copy[:, pivot_index]
        while_cont = 1
        cum_hist = np.zeros(M_input.shape[0])

        while while_cont == 1:
            while_cont = 0
            for j in Sremaining:
                cap_i = pivot_vector * M_copy[:, j]
                min_vec_size = np.min([np.sum(M_copy[:, j]), np.sum(pivot_vector)])
                cap_size = np.sum(cap_i)

                if cap_size / min_vec_size > overlap_coeff:
                    cum_hist = cum_hist + M_copy[:, j].astype(int)
                    while_cont = 1
                    Sremaining.remove(j)

        cnumT = np.floor(cum_hist.max() / hist_coeff)
        cum_vector = cum_hist > cnumT

        for j in ISet:
            capj = cum_vector * M_copy[:, j]
            difj = M_copy[:, j] > capj
            if np.sum(capj) > np.sum(difj):
                M_copy[:, j] = capj

            else:
                M_copy[:, j] = difj

        M_copy[:, pivot_index] = cum_vector
        ISet.remove(pivot_index)

    bret = np.argwhere(M_copy.astype(int) > M_input.astype(int))
    pret = np.argwhere(np.argwhere(M_copy.astype(int) < M_input.astype(int)))
    return [bret, pret, M_copy]


def ReadFfile(df):
    M = df.values.astype(bool)
    return M


def ReadFileNA(df):
    M = df.values
    NA_position = np.argwhere(M == 3)
    for j in NA_position:
        M[j[0], j[1]] = 0
    return M.astype(bool)


def Estimated_Matrix(df):
    M = df.values.astype(float)
    for i in range(M.shape[1]):
        if np.sum(M[:, i] != 3) == 0:
            one_ratio = 0
        else:
            one_ratio = np.sum(M[:, i] == 1) / np.sum(M[:, i] != 3)
        for j in range(M.shape[0]):
            if M[j, i] == 3:
                M[j, i] = one_ratio
    return M


def ReadFasis(df):
    M = df.values
    return M


def postprocess_col(df_input, matrix_recons, pfn, pfp):
    M_noisy = ReadFasis(df_input.copy())
    M_nds = matrix_recons.copy().astype(bool)
    Mtemp = c_m_col(ReadFasis(df_input.copy()), M_nds, pc_fn=pfn, pc_fp=pfp)
    Mtemp2 = Mtemp.copy()
    d10min = np.sum(Mtemp < (M_noisy == 1))
    imp = 1
    while imp:
        Mtemp2 = c_m_row(ReadFasis(df_input.copy()), Mtemp2, pc_fn=pfn, pc_fp=pfp)
        Mtemp2 = c_m_col(ReadFasis(df_input.copy()), Mtemp2, pc_fn=pfn, pc_fp=pfp)

        d10c = np.sum(Mtemp2 < (M_noisy == 1))
        if d10c < d10min:
            d10min = d10c
            Mtemp = Mtemp2.copy()
        else:
            imp = 0
    return Mtemp


def f_d_col(node_piv, M_samples, p_fp=0.005, p_fn=0.1):
    M_nodes = M_samples.copy()
    D10 = (M_nodes.T == 0).astype(int).dot(node_piv == 1)
    D11 = (M_nodes.T == 1).astype(int).dot(node_piv == 1)
    D00 = (M_nodes.T == 0).astype(int).dot(node_piv == 0)
    D01 = (M_nodes.T == 1).astype(int).dot(node_piv == 0)
    distances = -np.multiply(np.power(p_fp, D10), np.power(1 - p_fp, D00))
    distances = np.multiply(distances, np.power(p_fn, D01))
    distances = np.multiply(distances, np.power(1 - p_fn, D11))
    return distances


def c_m_col(M_input, M_nodes, pc_fp=0.0001, pc_fn=0.1):
    M_out = M_input.copy()
    for i in range(M_input.shape[1]):
        pivot_v = M_input[:, i]
        distance_i = f_d_col(pivot_v, M_nodes, p_fn=pc_fn, p_fp=pc_fp)
        min_index = np.argmin(distance_i)
        M_out[:, i] = M_nodes[:, min_index]
    return M_out


def f_d_row(node_piv, M_samples, p_fp=0.005, p_fn=0.1):
    M_nodes = M_samples.copy()
    D10 = (M_nodes == 0).astype(int).dot(node_piv == 1)
    D11 = (M_nodes == 1).astype(int).dot(node_piv == 1)
    D00 = (M_nodes == 0).astype(int).dot(node_piv == 0)
    D01 = (M_nodes == 1).astype(int).dot(node_piv == 0)
    distances = -np.multiply(np.power(p_fp, D10), np.power(1 - p_fp, D00))
    distances = np.multiply(distances, np.power(p_fn, D01))
    distances = np.multiply(distances, np.power(1 - p_fn, D11))
    return distances


def c_m_row(M_input, M_nodes, pc_fp=0.0001, pc_fn=0.1):
    M_out = M_input.copy()
    for i in range(M_input.shape[0]):
        pivot_v = M_input[i, :]
        distance_i = f_d_row(pivot_v, M_nodes, p_fn=pc_fn, p_fp=pc_fp)
        min_index = np.argmin(distance_i)
        M_out[i, :] = M_nodes[min_index, :]
    return M_out
