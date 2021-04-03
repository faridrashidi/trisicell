# -*- coding: utf-8 -*-
#!/usr/bin/env python
# -*- coding: utf-8 -*-


# =========================================================================================
# Written by : Can Kizilkale (CanKizilkale@lbl.gov) and Aydin Buluc (abuluc@lbl.gov)
# Last Update: Jan 31, 2021
# =========================================================================================

"""
Created on Thu Feb 13 19:18:26 2020

@author: Can KIZILKALE
"""

# Approximate Error Correction for phylogenetic tree matrices.
# Can Kizilkale

import argparse
import copy
import itertools
import multiprocessing as mp
import os
import time
from argparse import ArgumentParser
from datetime import datetime
from multiprocessing import Process, Queue

import networkx as nx
import numpy as np
import pandas as pd

os.environ["PATH"] += os.pathsep + "C:/Program Files (x86)/Graphviz2.38/bin/"
np.set_printoptions(threshold=np.inf)

# ---------------------READ ME ---------------------------
# Call this function like:
#
# \python huntress.py "Noisy_matrix_filename" "Output_filename" --nofcpus 8 --algorithmchoice "FPNA" --fn_fpratio 100 --fp_coeff 0.0001 --fn_coeff 0.1
# Output is written on Output_filename.CFMATRIX
#
# The path and name of the the noisy matrix is given in Noisy_matrix_filename
# The reconstructed error free matrix is written in Output_filename with the extension ".CFMatrix"
# The optional inputs are as follows:
# --nofcpus defines the number of cpus to be used for tuning in parallel. Default is 7
# --algorithmchoice defines the version of the algorithm to be used.
#           = "FN" for matrices that only have false negatives
#           = "FPNA" for matrices that have false positives , false negatives and NA (entries that could not be read) entries. These entries must be given as 3 in the input matrix
# --fn_fpratio is the ratio of the weights of 0->1 switches over 1->0 switches that is used by the algorithm to tune the parameters.
# Default=100
# --fp_coeff false positive probability coefficient used for postprocessing.
# Default: 0.0001
# --fn_coeff false negative probability coefficient used for postprocessing.
# Default: 0.1
#
# ------------------- END README ---------------------------------------

prfp = 5


def Reconstruct(
    input_file,
    output_file,
    Algchoice="FPNA",
    auto_tune=1,
    overlapp_coeff=0.3,
    hist_coeff=20,
    postprocessing=0,
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
    Flog = open(output_file + ".LOG", "a+")
    matrix_input = ReadFileNA(input_file)
    #    matrix_input=preproc_row(matrix_input,0.5)
    matrix_input_raw = ReadFasis(input_file)
    matrix_NA_est = Estimated_Matrix(input_file)
    print(np.sum(matrix_input), np.sum(matrix_NA_est), file=Flog)
    h_range = [x for x in range(1, 100, 2)]
    oc_range = [x / 10 for x in range(1, 5)]
    tune_var = [x for x in itertools.product(h_range, oc_range)]

    running_time = 0

    if Algchoice == "FN":
        s_time = time.time()
        matrix_recons = greedyPtreeNew(matrix_input.astype(bool))[1]
        e_time = time.time()
        running_time = e_time - s_time
        WriteTfile(output_file, matrix_recons, input_file)

    if Algchoice == "FPNA" and auto_tune == 0:

        s_time = time.time()
        apprx_ordr = sum(matrix_NA_est)
        matrix_recons = greedyPtreeNA(
            matrix_input.astype(bool), apprx_ordr, overlapp_coeff, hist_coeff
        )[2]
        e_time = time.time()
        running_time = e_time - s_time
        output_file = output_file + "_optH_{}TEMP.CFMatrix".format(hist_coeff)
        WriteTfile(output_file, matrix_recons, input_file)
    #        print(" difference ",np.sum(matrix_recons!=matrix_input))

    if Algchoice == "FPNA" and auto_tune == 1:
        proc_size = np.ceil(len(tune_var) / n_proc).astype(int)
        print(proc_size)
        cpu_range = []
        #        q=Queue()

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
                    i,
                    input_file,
                    output_file,
                ),
            )
            for i in range(len(cpu_range))
        ]
        for i in p:
            i.start()
        for i in p:
            i.join()
        ret = []
        while not q.empty():
            print(" Reading the queue ...")
            ret.append(q.get())

        [m_r, d_min] = ret[0]
        matrix_recons = ReadFfile(m_r)
        for i in range(1, len(ret)):
            [m, d] = ret[i]
            print(d_min)
            if d < d_min:
                matrix_recons = ReadFfile(m)
                print(d_min)
                d_min = d

        print("closed", q.empty())

        e_time = time.time()
        running_time = e_time - s_time
        print(running_time)
        output_file = output_file + "_optH_TEMP.CFMatrix"
        WriteTfile(output_file, matrix_recons, input_file)
        print(
            " 1->0 : ",
            np.sum(matrix_recons < matrix_input),
            " 0->1 : ",
            np.sum(matrix_recons > matrix_input_raw),
            " NA->1 : ",
            np.sum(matrix_recons > matrix_input)
            - np.sum(matrix_recons > matrix_input_raw),
        )

    Flog.close()
    postprocess_col(
        input_file, output_file, pfn=post_fn, pfp=post_fp
    )  # Column based post processing, seem to give the best result.

    return running_time


def Auto_fnfp(
    q, tune_ran, m_input, m_NA_est, m_input_raw, fnfp, fnc, procid, in_file, out_file
):  # multiprocessing worker unit, takes a tuning range and tries everything in between
    apprx_ordr = sum(m_NA_est)

    print("An instance of auto tune has started...", procid)
    matrix_recon = greedyPtreeNA(
        m_input.astype(bool), apprx_ordr, tune_ran[0][1], tune_ran[0][0]
    )[2]
    matrix_rec_Temp = deleteNas(m_input_raw, matrix_recon)
    n_10 = np.sum(matrix_rec_Temp < m_input)
    n_01 = np.sum(matrix_rec_Temp > m_input)

    #        distance=np.square(n_10)+n_01
    distance = fnfp * n_10 + fnc * n_01

    h_current = tune_ran[0][0]
    oc_current = tune_ran[0][1]
    for [h_i, overlapp_coeff] in tune_ran:
        matrix_rec_i = greedyPtreeNA(
            m_input.astype(bool), apprx_ordr, overlapp_coeff, h_i
        )[2]
        matrix_rec_Temp = deleteNas(m_input_raw, matrix_rec_i)
        n_10 = np.sum(matrix_rec_Temp < m_input, dtype="int64")
        n_01 = np.sum(matrix_rec_Temp > m_input, dtype="int64")
        #        print("Opt h = ",h_current,"Opt_O = ",oc_current,"H,OC coefficient = ",h_i,overlapp_coeff," 01 switches : ",n_01,"10 switches = ", n_10 , "best :",distance)
        print(procid, h_i, overlapp_coeff)
        #           distance_i=np.square(n_10)+n_01
        distance_i = fnfp * n_10 + fnc * n_01
        #        distance_i=-fnfp**n_10*fnc**n_01
        #            if n_10<np.sum(matrix_recons<matrix_input):
        if distance_i < distance:
            matrix_recon = matrix_rec_i.copy()
            distance = distance_i
            #                WriteTfile(output_file,matrix_recons,input_file)
            h_current = h_i
            oc_current = overlapp_coeff
    print(procid, "th process finished", q.full(), h_current, oc_current)
    output_file = out_file + "_TEMP_{}.CFMatrix".format(procid)
    WriteTfile(output_file, matrix_recon, in_file)
    q.put([output_file, distance])


def deleteNas(M_in, M_out):
    M_o = M_out.copy()
    NA_position = np.argwhere(M_in == 3)
    # print("number of NA : ",len(NA_position))
    for j in NA_position:
        # print(M[j[0],j[1]])
        M_o[j[0], j[1]] = 0

    return M_o


def deletemutations(
    M_in, M_out
):  # Finds the rows which are too far from the input, just replicates them to their closest neighbour.
    x = M_in.shape
    M_return = M_out.copy()
    treshold_cut = 0.5
    dt = np.divide(sum(M_in), sum(M_out)) <= treshold_cut
    for i in range(x[1]):
        if dt[i] == 1:
            M_return[:, i] = np.zeros(x[0])

    return M_return


def precombination(M_in):
    x = M_in.shape
    trsh = 0.22
    M_return = M_in.copy()
    for i in range(x[0]):
        for j in range(i + 1, x[0]):
            rij = np.sum(M_in[i, :] != M_in[j, :]) / np.sum(M_in[i, :] + M_in[j, :])
            print(rij)
            if rij < trsh:
                cupi = M_in[i, :] + M_in[j, :]
                M_return[i, :] = cupi
                M_return[j, :] = cupi
                print("combined ", i, j)
    return M_return


def findclosestinter(i, M_input):
    vec_int = np.zeros(M_input.shape[1]).astype(bool)
    for j in range(M_input.shape[0]):
        if (
            i != j
            and np.sum(M_input[j, :] > M[i, :]) > 0
            and np.sum(M_input[j, :] * M[i, :]) > 0
        ):
            if np.sum(vec_int):
                vec_int = M_input[j, :]
            else:
                vec_int = vec_int * M_input[j, :]
    return vec_int


# Both algorithms take inut matrices of BOOL type.
def greedyPtreeNew(
    M_input,
):  # very greedy algorithm that constructs a ptree matrix from M_inputs by adding 1's
    # M_input has to be of type BOOLEAN !
    # Returns
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    #    M_apprx=greedy_row_rec(M_input)[1]
    M_copy = M_input.copy()

    ISet1 = np.argsort(sum(M_copy))
    #    ISet1=np.argsort(sum(M_apprx))
    #    ISet1=np.argsort(sum(M_rowdene))
    ISet = []
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret = []  # Location of detected false negatives
    print(M_copy.shape, len(ISet))
    while len(ISet) > 1:
        # print("pivoting column", i,bret)
        pivot_index = ISet[-1]  # index of the pivot vector
        Sremaining = ISet.copy()
        pivot_vector = M_copy[
            :, pivot_index
        ]  # vector used for pivoting the current iteration
        cum_vector = np.copy(pivot_vector)  # holds the union vector
        while_cont = 1

        while while_cont == 1:  # main loop
            while_cont = 0

            for j in Sremaining:
                cap_j = (
                    cum_vector * M_copy[:, j]
                )  # intersection of the pivot and the jth column

                if np.any(
                    cap_j
                ):  # continue as long as there is a column having non-empty intersection
                    cum_vector = cum_vector + M_copy[:, j]
                    while_cont = 1
                    Sremaining.remove(j)
        #                    print("pivot ", pivot_index, j, sum(cap_j),sum(cum_vector),sum(M_copy[:,j]))
        #        print(len(Si))
        #        print(sum(pivot_vector),sum(cum_vector))
        M_copy[:, pivot_index] = cum_vector
        ISet.remove(pivot_index)
        bret = np.argwhere(M_copy.astype(int) > M_input.astype(int))
    return [bret, M_copy]


def greedyPtreeNA(
    M_input, approx_order, oc, hc
):  # Modified Greedy algorithm for NA and false positive values.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff = oc  # 0.1
    hist_coeff = hc  # 25

    M_copy = M_input.copy()
    #    ISet1=np.argsort(sum((M_copy.T.dot(M_copy)>0).T))    # Set of indices of columns sorted wrt. number of overlapping vectors
    #    ISet1=np.argsort(sum(M_copy))         # Set of indices of columns sorted wrt. the number of ones they contain
    # ISet1=np.argsort(sum(M_NA))
    ISet1 = np.argsort(approx_order)
    ISet = []
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret = []  # Location of detected false negatives
    pret = []  # Location of detected false positives

    while len(ISet) > 1:
        #        pivot_index=ISet[np.argmax(sum((M_copy[:,ISet].T.dot(M_copy[:,ISet])>0).T))]
        pivot_index = ISet[-1]  # index of the pivot vector
        Sremaining = ISet.copy()  # set of indices that are not included in the union
        pivot_vector = M_copy[
            :, pivot_index
        ]  # vector used for pivoting the current iteration
        cum_vector = np.copy(pivot_vector)  # holds the union vector
        while_cont = 1
        cum_hist = np.zeros(M_input.shape[0])  # holds the histogram for the union

        while while_cont == 1:  # Continue uniting vectors until no candidate remains
            while_cont = 0
            for j in Sremaining:
                cap_i = pivot_vector * M_copy[:, j]
                min_vec_size = np.min([np.sum(M_copy[:, j]), np.sum(pivot_vector)])
                cap_size = np.sum(cap_i)
                #               print(np.sum(M_copy[:,j]),np.sum(pivot_vector),np.sum(cap_i),cap_size, min_vec_size,np.floor(cap_size/min_vec_size))

                if (
                    cap_size / min_vec_size > overlap_coeff
                ):  # we check if the columns have a meaningful overlap
                    cum_hist = cum_hist + M_copy[:, j].astype(int)
                    while_cont = 1  # found an overlapping vector so we keep on going
                    Sremaining.remove(j)  # united vector is removed from the search set

        cnumT = np.floor(
            cum_hist.max() / hist_coeff
        )  # the elements that repeated few times are considered to be false positives
        cum_vector = cum_hist > cnumT

        pivot_est = pivot_index

        ncap = np.sum(cum_vector * M_copy[:, pivot_est])
        #        ndel=np.sum(M_copy[:,j]!=capj)
        for j in ISet:  # clean up the false positives wrt. the established pivot
            capj = cum_vector * M_copy[:, j]  # intersection of union with column j
            #            difj=M_copy[:,j]!=capj                        # difference of column j from the union
            difj = M_copy[:, j] > capj  # difference of column j from the union
            if np.sum(capj) > np.sum(difj):
                M_copy[:, j] = capj

            else:
                M_copy[:, j] = difj

        M_copy[
            :, pivot_index
        ] = cum_vector  # correcting the false negatives in the pivot
        ISet.remove(pivot_index)  # removing the pivot from the search space
    #        print(len(ISet),np.sum(M_copy[:,pivot_index]),np.sum(M_input[:,pivot_index]))

    bret = np.argwhere(M_copy.astype(int) > M_input.astype(int))
    pret = np.argwhere(np.argwhere(M_copy.astype(int) < M_input.astype(int)))
    return [bret, pret, M_copy]


def ReadFfile(filename):  # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep="\t", index_col=0)
    M = df.values.astype(bool)
    return M


def ReadFile(filename):  # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep="\t", index_col=0)
    M = df.values
    return M


def ReadFileNA(filename):  # reads the file and fills the NA with 0's.
    df = pd.read_csv(filename, sep="\t", index_col=0)
    M = df.values

    NA_position = np.argwhere(M == 3)
    # print("number of NA : ",len(NA_position))
    for j in NA_position:
        # print(M[j[0],j[1]])
        M[j[0], j[1]] = 0

    # print(sum(sum(M)))
    return M.astype(bool)


def Estimated_Matrix(
    filename,
):  # Creates an estimate of the matrix such that each element is given the expectation wrt the column 1/0 frequencies.
    df = pd.read_csv(filename, sep="\t", index_col=0)
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


def WriteTfile(
    filename, matrix, filename2
):  # writes matrix output as an integer matrix
    df_input = pd.read_csv(filename2, sep="\t", index_col=0)
    matrix_output = matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"
    df_output.to_csv(filename, sep="\t")


##########################################################################3
############# TESTING FUNCTIONS ##########################################
def isPtree(matrix_in):  # brute force check if matrix_in is a pTree
    M = matrix_in.astype(bool)
    for j in range(M.shape[1]):
        for i in range(j, M.shape[1]):
            cap = M[:, i] * M[:, j]
            cap_size = np.sum(cap)
            Mi_size = np.sum(M[:, i])
            Mj_size = np.sum(M[:, j])
            if cap_size != 0:
                if cap_size != Mi_size:
                    if cap_size != Mj_size:
                        return False

    print("Seems to be a PTree ...")
    return True


def compareAD(M1, M2):  # M1 is the ground truth
    error_pairs = []
    n_adpairs = 0
    for i in range(M1.shape[1]):
        #        print(i)
        for j in range(i, M1.shape[1]):
            cap1 = M1[:, i] * M1[:, j]
            cap2 = M2[:, i] * M2[:, j]
            if np.sum(cap1) > 0 and np.sum(M1[:, i]) != np.sum(M1[:, j]):
                n_adpairs = n_adpairs + 1
                if np.sum(cap2) == 0:
                    error_pairs.append([i, j])
                else:
                    if np.sum(M1[:, j]) > np.sum(M1[:, i]) and np.sum(
                        M2[:, j]
                    ) <= np.sum(M2[:, i]):
                        error_pairs.append([i, j])
                    else:
                        if np.sum(M1[:, i]) > np.sum(M1[:, j]) and np.sum(
                            M2[:, i]
                        ) <= np.sum(M2[:, j]):
                            error_pairs.append([i, j])
                        # print(i,j,sum(M1[:,i]),sum(M1[:,j]),sum(M2[:,i]),sum(M2[:,j]))
    print(
        "Number of AD pairs = ",
        n_adpairs,
        "errors : ",
        len(error_pairs),
        "AD score = ",
        1 - len(error_pairs) / n_adpairs,
    )
    return error_pairs


def compareDF(M_orj, M_rec):
    error_pairs = []
    d_pairs = 0
    for i in range(M_orj.shape[1]):
        for j in range(i, M_orj.shape[1]):
            cap1 = M_orj[:, i] * M_orj[:, j]
            cap2 = M_rec[:, i] * M_rec[:, j]
            if np.sum(cap1) == 0:
                d_pairs = d_pairs + 1
                if np.sum(cap2) > 0:
                    error_pairs.append([i, j])
    print(
        "Number of Diff pairs = ",
        d_pairs,
        "errors :",
        len(error_pairs),
        "score :",
        1 - len(error_pairs) / d_pairs,
    )
    return


def ReadFasis(filename):  # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep="\t", index_col=0)
    M = df.values
    return M


def compute_fnfp(M_n, M_r):
    n_01 = 0
    n_10 = 0
    n_1 = 0
    n_0 = 0
    for x in np.argwhere(M_n < 3):
        if M_r[x[0], x[1]] == 0:
            n_0 = n_0 + 1
        if M_r[x[0], x[1]] == 1:
            n_1 = n_1 + 1
        if M_n[x[0], x[1]] > M_r[x[0], x[1]]:

            n_10 = n_10 + 1

        if M_n[x[0], x[1]] < M_r[x[0], x[1]]:
            n_01 = n_01 + 1
            n_1 = n_1 + 1
    print("computed fn :", n_01 / n_1, " fp : ", n_10 / n_0)
    return [n_01 / n_1, n_10 / n_0]


def find_dist(node_piv, M_samples):
    #    distances=np.zeros(M_nodes.shape[0])
    #    for i in range(M_nodes.shape[0]):
    ##        distances[i]=np.sum(M_nodes[i,:]!=node_piv)
    ##        distances[i]=np.sum(M_nodes[i,:]<node_piv)
    ##        distances[i]=np.linalg.norm(M_nodes[i,:]-node_piv)
    #        distances[i]=0
    #        d_10=0
    #        d_01=0
    #        for j in range(M_nodes.shape[1]):
    #
    #            if node_piv[j]!=3:
    #                if node_piv[j]>M_nodes[i,j]:
    #                    d_10=d_10+1
    #                if node_piv[j]<M_nodes[i,j]:
    #                    d_01=d_01+1
    #        distances[i]= np.square(d_10) +d_01
    ##            else:
    ##                distances[i]=1-M_nodes[i,j]+distances[i]
    M_nodes = M_samples.copy()
    for j in range(M_nodes.shape[1]):
        if node_piv[j] == 3:
            node_piv[j] = 0
            M_nodes[:, j] = 0 * M_nodes[:, j]
    distances = np.zeros(M_nodes.shape[0])
    for i in range(M_nodes.shape[0]):
        d_10 = np.sum(node_piv > M_nodes[i, :], dtype="int64")
        d_01 = np.sum(node_piv < M_nodes[i, :], dtype="int64")
        distances[i] = np.square(d_10) + d_01

    return distances


def find_dist_col(node_piv, M_samples):
    M_nodes = M_samples.copy()
    for j in range(M_nodes.shape[0]):
        if node_piv[j] == 3:
            node_piv[j] = 0
            M_nodes[j, :] = 0 * M_nodes[j, :]
    distances = np.zeros(M_nodes.shape[1])
    for i in range(M_nodes.shape[1]):
        d_10 = np.sum(node_piv > M_nodes[:, i], dtype="int64")
        d_01 = np.sum(node_piv < M_nodes[:, i], dtype="int64")
        #        distances[i]=np.square(d_10) + d_01
        #        distances[i]=d_10 + d_01
        distances[i] = prfp * d_10 + d_01
        #        distances[i]=np.sum(M_nodes[:,i]*node_piv)/np.sqrt(np.sum(M_nodes[:,i])*np.sum(node_piv))
        distances[i] = -((0.005) ** d_10) * (0.1) ** d_01
    return distances


def closest_matrix(M_input, M_nodes, M_rec):
    M_out = M_input.copy()
    for i in range(M_input.shape[0]):
        pivot_v = M_input[i, :]
        distance_i = find_dist(pivot_v, M_nodes)

        min_index = np.argmin(distance_i)
        #        print("Current distance ",np.sum(M_input[i,:]>M_rec[i,:]),"minimum possible ",np.min(distance_i),min_index)
        M_out[i, :] = M_nodes[min_index, :]
    #        if np.sum(M_input[i,:]>M_rec[i,:])>2*np.min(distance_i):
    #            min_index=np.argmin(distance_i)
    #            M_out[i,:]=M_nodes[min_index,:]
    #        else:
    #            M_out[i,:]=M_rec[i,:]
    #        print("M input :", sum(M_input[i,:]),sum(M_out[i,:]))
    return M_out


def closest_matrix_col(M_input, M_nodes, M_rec):
    M_out = M_input.copy()
    for i in range(M_input.shape[1]):
        pivot_v = M_input[:, i]
        distance_i = find_dist_col(pivot_v, M_nodes)

        min_index = np.argmin(distance_i)
        #        print("Old New difference ",np.sum(M_nodes[:,i]!=M_nodes[:,min_index]),i)
        M_out[:, i] = M_nodes[:, min_index]

    return M_out


def postprocess_col(input_file, out_file, pfn, pfp):
    s = time.time()

    M_noisy = ReadFasis(input_file)
    M_n_copy = M_noisy.copy()
    M_nds = ReadFfile(out_file)

    Mtemp = c_m_col(ReadFasis(input_file), M_nds, pc_fn=pfn, pc_fp=pfp)
    Mtemp2 = Mtemp.copy()
    #    Mtemp=c_m_row(ReadFasis(input_file),Mtemp,1)
    d10min = np.sum(Mtemp < (M_noisy == 1))
    d10c = d10min
    imp = 1
    while imp:
        Mtemp2 = c_m_row(ReadFasis(input_file), Mtemp2, pc_fn=pfn, pc_fp=pfp)
        Mtemp2 = c_m_col(ReadFasis(input_file), Mtemp2, pc_fn=pfn, pc_fp=pfp)

        #        Mtemp=c_m_col(ReadFasis(input_file),Mtemp,1)
        d10c = np.sum(Mtemp2 < (M_noisy == 1))
        print(d10c)
        if d10c < d10min:
            d10min = d10c
            Mtemp = Mtemp2.copy()
        else:
            imp = 0

    #    M_postprocessed=closest_matrix_col(M_noisy,M_nds,ReadFfile(out_file))
    M_postprocessed = Mtemp
    processed_file = out_file[:-13] + ".CFMatrix"
    #    print("Writing to file ",processed_file,file=Flog)
    WriteTfile(processed_file, M_postprocessed, input_file)
    #    print(np.sum(M_noisy),np.sum(M_noisy==1))
    e = time.time()
    print(
        "Postprocessed 1->0 : ",
        np.sum(M_postprocessed < (M_n_copy == 1)),
        " 0->1 : ",
        np.sum(M_postprocessed > M_n_copy),
        " NA->1 : ",
        np.sum(M_postprocessed > (M_n_copy == 1)) - np.sum(M_postprocessed > M_n_copy),
    )
    print("Post processing time : ", e - s)


#    draw_tree(processed_file)


def preproc_row(M_o, c=0.8):
    M_pre = M_o.copy()
    for i in range(M_o.shape[0]):
        for j in range(i + 1, M_o.shape[0]):
            prdct = np.sum(M_o[i, :] * M_o[j, :]) / np.sqrt(
                np.sum(M_o[i, :]) * np.sum(M_o[j, :])
            )
            if prdct > c:
                print(i, j)
                M_pre[i, :] = M_o[i, :] + M_o[j, :]
                M_pre[j, :] = M_o[i, :] + M_o[j, :]
    return M_pre


def f_d_col(node_piv, M_samples, p_fp=0.005, p_fn=0.1):
    M_nodes = M_samples.copy()
    D10 = (M_nodes.T == 0).astype(int).dot(node_piv == 1)
    D11 = (M_nodes.T == 1).astype(int).dot(node_piv == 1)
    D00 = (M_nodes.T == 0).astype(int).dot(node_piv == 0)
    D01 = (M_nodes.T == 1).astype(int).dot(node_piv == 0)
    distances = np.zeros(M_nodes.shape[1])
    distances = -np.multiply(np.power(p_fp, D10), np.power(1 - p_fp, D00))
    distances = np.multiply(distances, np.power(p_fn, D01))
    distances = np.multiply(distances, np.power(1 - p_fn, D11))
    #    for i in range(M_nodes.shape[1]):
    #        d_10=D10[i]
    #        d_11=D11[i]
    #        d_00=D00[i]
    #        d_01=D01[i]
    ##        d_10=np.sum((M_nodes[:,i]==0)*(node_piv==1),dtype='int64')
    ##        d_11=np.sum((M_nodes[:,i]==1)*(node_piv==1),dtype='int64')
    ##        d_00=np.sum((M_nodes[:,i]==0)*(node_piv==0),dtype='int64')
    ##        d_01=np.sum((M_nodes[:,i]==1)*(node_piv==0),dtype='int64')
    #
    #        distances[i]=-(p_fp)**d_10*(1-p_fp)**d_00*(p_fn)**d_01*(1-p_fn)**d_11
    return distances


def c_m_col(M_input, M_nodes, pc_fp=0.0001, pc_fn=0.1):
    M_out = M_input.copy()

    for i in range(M_input.shape[1]):
        pivot_v = M_input[:, i]

        distance_i = f_d_col(pivot_v, M_nodes, p_fn=pc_fn, p_fp=pc_fp)

        min_index = np.argmin(distance_i)
        #        print("Old New difference ",np.sum(M_nodes[:,i]!=M_nodes[:,min_index]),i)
        M_out[:, i] = M_nodes[:, min_index]

    return M_out


def f_d_row(node_piv, M_samples, p_fp=0.005, p_fn=0.1):
    M_nodes = M_samples.copy()

    distances = np.zeros(M_nodes.shape[0])
    D10 = (M_nodes == 0).astype(int).dot(node_piv == 1)
    D11 = (M_nodes == 1).astype(int).dot(node_piv == 1)
    D00 = (M_nodes == 0).astype(int).dot(node_piv == 0)
    D01 = (M_nodes == 1).astype(int).dot(node_piv == 0)

    distances = -np.multiply(np.power(p_fp, D10), np.power(1 - p_fp, D00))
    distances = np.multiply(distances, np.power(p_fn, D01))
    distances = np.multiply(distances, np.power(1 - p_fn, D11))
    #    for i in range(M_nodes.shape[0]):
    #        d_10=D10[i]
    #        d_11=D11[i]
    #        d_00=D00[i]
    #        d_01=D01[i]
    ##
    ##        d_10=np.sum((M_nodes[i,:]==0)*(node_piv==1),dtype='int64')
    ##        d_11=np.sum((M_nodes[i,:]==1)*(node_piv==1),dtype='int64')
    ##        d_00=np.sum((M_nodes[i,:]==0)*(node_piv==0),dtype='int64')
    ##        d_01=np.sum((M_nodes[i,:]==1)*(node_piv==0),dtype='int64')
    #
    #        distances[i]=-(p_fp)**d_10*(1-p_fp)**d_00*(p_fn)**d_01*(1-p_fn)**d_11
    return distances


def c_m_row(M_input, M_nodes, pc_fp=0.0001, pc_fn=0.1):
    M_out = M_input.copy()

    for i in range(M_input.shape[0]):
        pivot_v = M_input[i, :]
        distance_i = f_d_row(pivot_v, M_nodes, p_fn=pc_fn, p_fp=pc_fp)

        min_index = np.argmin(distance_i)
        #        print("Old New difference ",np.sum(M_nodes[:,i]!=M_nodes[:,min_index]),i)
        M_out[i, :] = M_nodes[min_index, :]

    return M_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile")
    parser.add_argument("outputfile")
    parser.add_argument("--nofcpus", default=7, type=int, nargs="?")
    parser.add_argument("--algorithmchoice", default="FPNA", nargs="?")
    parser.add_argument("--fn_fpratio", default=51, type=int, nargs="?")
    parser.add_argument("--fp_coeff", default=0.00001, type=float, nargs="?")
    parser.add_argument("--fn_coeff", default=0.1, type=float, nargs="?")
    args = parser.parse_args()

    fn_conorm = 0.1
    fp_conorm = fn_conorm * args.fp_coeff / args.fn_coeff
    #    fnfp_conorm = fn_conorm / fp_conorm
    #    Reconstruct(args.inputfile,args.outputfile,Algchoice=args.algorithmchoice,n_proc=args.nofcpus,fnfp=args.fn_fpratio,post_fn=args.fn_coeff,post_fp=args.fp_coeff)
    Reconstruct(
        args.inputfile,
        args.outputfile,
        Algchoice=args.algorithmchoice,
        n_proc=args.nofcpus,
        fnfp=args.fn_fpratio,
        post_fn=fn_conorm,
        post_fp=fp_conorm,
    )
