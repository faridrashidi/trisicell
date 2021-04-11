import contextlib
import datetime
import functools
import multiprocessing
import os
import shutil
import tempfile
import time

import joblib
import numpy as np
import pandas as pd
import pkg_resources

import trisicell as tsc


def log_input(df_in):
    size = df_in.shape[0] * df_in.shape[1]
    tsc.logg.info(f"input -- size: {df_in.shape[0]}x{df_in.shape[1]}")
    tsc.logg.info(
        f"input -- 0: {np.sum(df_in.values == 0)}#,"
        f" {100*np.sum(df_in.values == 0)/size:.1f}%"
    )
    tsc.logg.info(
        f"input -- 1: {np.sum(df_in.values == 1)}#,"
        f" {100*np.sum(df_in.values == 1)/size:.1f}%"
    )
    tsc.logg.info(
        f"input -- NA: {np.sum(df_in.values == 3)}#,"
        f" {100*np.sum(df_in.values == 3)/size:.1f}%"
    )
    tsc.logg.info(f"input -- CF: {is_conflict_free_gusfield(df_in)}")


def log_output(df_out, running_time):
    size = df_out.shape[0] * df_out.shape[1]
    tsc.logg.info(f"output -- size: {df_out.shape[0]}x{df_out.shape[1]}")
    tsc.logg.info(
        f"output -- 0: {np.sum(df_out.values == 0)}#,"
        f" {100*np.sum(df_out.values == 0)/size:.1f}%"
    )
    tsc.logg.info(
        f"output -- 1: {np.sum(df_out.values == 1)}#,"
        f" {100*np.sum(df_out.values == 1)/size:.1f}%"
    )
    tsc.logg.info(
        f"output -- NA: {np.sum(df_out.values == 3)}#,"
        f" {100*np.sum(df_out.values == 3)/size:.1f}%"
    )
    tsc.logg.info(f"output -- CF: {is_conflict_free_gusfield(df_out)}")
    tsc.logg.info(
        f"output -- time: {running_time:.1f}s"
        f" ({datetime.timedelta(seconds=running_time)})"
    )


def log_flip(df_in, df_out):
    flips_0_1, flips_1_0, flips_na_0, flips_na_1 = count_flips(
        df_in.values, df_out.values, 3
    )
    fn_rate, fp_rate, na_rate = infer_rates(df_in.values, df_out.values, 3)
    tsc.logg.info(f"flips -- #0->1: {flips_0_1}")
    tsc.logg.info(f"flips -- #1->0: {flips_1_0}")
    tsc.logg.info(f"flips -- #NA->0: {flips_na_0}")
    tsc.logg.info(f"flips -- #NA->1: {flips_na_1}")
    tsc.logg.info(f"rates -- FN: {fn_rate:.3f}")
    tsc.logg.info(f"rates -- FP: {fp_rate:.8f}")
    tsc.logg.info(f"rates -- NA: {na_rate:.3f}")


def calc_score_tree(df_in, df_out, alpha, beta):
    if alpha == 0 or beta == 0:
        return None
    columns = np.intersect1d(df_in.columns, df_out.columns)
    indices = np.intersect1d(df_in.index, df_out.index)
    D = df_in.loc[indices, columns].values
    E = df_out.loc[indices, columns].values
    removedMutations = []
    objective = 0
    for j in range(len(columns)):
        numZeros = 0
        numOnes = 0
        for i in range(len(indices)):
            if D[i, j] == 0:
                numZeros += 1
                objective += np.log(beta / (1 - alpha)) * E[i, j]
            elif D[i, j] == 1:
                numOnes += 1
                objective += np.log((1 - beta) / alpha) * E[i, j]
        objective += numZeros * np.log(1 - alpha)
        objective += numOnes * np.log(alpha)
        if j in removedMutations:
            objective -= numZeros * np.log(1 - alpha) + numOnes * (
                np.log(alpha) + np.log((1 - beta) / alpha)
            )
    tsc.logg.info(f"score -- NLL: {-objective}")
    return None


def stat(df_in, df_out, alpha, beta, running_time):
    log_input(df_in)
    log_output(df_out, running_time)
    log_flip(df_in, df_out)
    calc_score_tree(df_in, df_out, alpha, beta)


def get_param(filename):
    data = {}
    # simNo_2-s_7-m_20-h_1-minVAF_0.1-ISAV_0-n_10-fp_0-fn_0.1-na_0-d_0-l_1000000.SC
    dirname, basename = dir_base(filename)
    data["simNo"] = int(basename.split("-")[0].split("_")[1])
    data["s"] = int(basename.split("-")[1].split("_")[1])
    data["m"] = int(basename.split("-")[2].split("_")[1])
    data["h"] = int(basename.split("-")[3].split("_")[1])
    data["minVAF"] = float(basename.split("-")[4].split("_")[1])
    data["ISAV"] = int(basename.split("-")[5].split("_")[1])
    data["n"] = int(basename.split("-")[6].split("_")[1])
    data["fp"] = float(basename.split("-")[7].split("_")[1])
    data["fn"] = float(basename.split("-")[8].split("_")[1])
    data["na"] = float(basename.split("-")[9].split("_")[1])
    data["d"] = int(basename.split("-")[10].split("_")[1])
    last = basename.split("-")[11]
    if "." in last:
        data["l"] = int(last.split(".")[0].split("_")[1])
    else:
        data["l"] = int(last.split("_")[1])
    return data


def count_flips(I, O, na_value=3):
    flips_0_1 = 0
    flips_1_0 = 0
    flips_na_0 = 0
    flips_na_1 = 0
    n, m = I.shape
    for i in range(n):
        for j in range(m):
            if I[i, j] == 0 and O[i, j] == 1:
                flips_0_1 += 1
            elif I[i, j] == 1 and O[i, j] == 0:
                flips_1_0 += 1
            elif I[i, j] == na_value and O[i, j] == 0:
                flips_na_0 += 1
            elif I[i, j] == na_value and O[i, j] == 1:
                flips_na_1 += 1
    return flips_0_1, flips_1_0, flips_na_0, flips_na_1


def infer_rates(I, O, na_value=3):
    flips_0_1, flips_1_0, flips_na_0, flips_na_1 = count_flips(I, O, na_value)
    fn_rate = flips_0_1 / ((O == 1) & (I != na_value)).sum()
    fp_rate = flips_1_0 / ((O == 0) & (I != na_value)).sum()
    na_rate = (flips_na_1 + flips_na_0) / I.size
    return fn_rate, fp_rate, na_rate


def is_conflict_free(D, na_value=3):
    conflict_free = True
    for p in range(D.shape[1]):
        for q in range(p + 1, D.shape[1]):
            oneone = False
            zeroone = False
            onezero = False
            for r in range(D.shape[0]):
                if D[r, p] == 1 and D[r, q] == 1:
                    oneone = True
                if D[r, p] == 0 and D[r, q] == 1:
                    zeroone = True
                if D[r, p] == 1 and D[r, q] == 0:
                    onezero = True
            if oneone and zeroone and onezero:
                conflict_free = False
                return conflict_free
    return conflict_free


def is_conflict_free_gusfield(df_in, na_value=3):
    """Check conflict-free criteria via Gusfield algorithm.

    This is an implementation of algorithm 1.1 in :cite:`Gusfield_1991`.

    The order of this algorithm is O(nm)
    where n is the number of cells and m is the number of mutations.

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.
    na_value : :obj:`int`, optional
        Missing value representation in the input data, by default 3.

    Returns
    -------
    :obj:`bool`
        A Boolean checking if the input conflict-free or not.

    Examples
    --------
    >>> sc = tsc.datasets.test()
    >>> tsc.ul.is_conflict_free_gusfield(sc)
    False

    See Also
    --------
    :func:`trisicell.ul.is_conflict_free`.
    """

    I = df_in.astype(int).values

    def sort_bin(a):
        b = np.transpose(a)
        b_view = np.ascontiguousarray(b).view(
            np.dtype((np.void, b.dtype.itemsize * b.shape[1]))
        )
        idx = np.argsort(b_view.ravel())[::-1]
        c = b[idx]
        return np.transpose(c), idx

    Ip = I.copy()
    # Ip[Ip == na_value] = 0
    O, idx = sort_bin(Ip)
    # tsc.logg.info(O, '\n')
    Lij = np.zeros(O.shape, dtype=int)
    for i in range(O.shape[0]):
        maxK = 0
        for j in range(O.shape[1]):
            if O[i, j] == 1:
                Lij[i, j] = maxK
                maxK = j + 1
    # tsc.logg.info(Lij, '\n')
    Lj = np.amax(Lij, axis=0)
    # tsc.logg.info(Lj, '\n')
    for i in range(O.shape[0]):
        for j in range(O.shape[1]):
            if O[i, j] == 1:
                if Lij[i, j] != Lj[j]:
                    return False  # , (idx[j], idx[Lj[j] - 1])
    return True  # , (None, None)


def tmpdir(prefix="trisicell.", suffix=".trisicell", dirname="."):
    return tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dirname)


def tmpfile(prefix="trisicell.", suffix=".trisicell", dirname="."):
    return tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dirname)[1]


def tmpdirsys(prefix="trisicell.", suffix=".trisicell", dirname="."):
    return tempfile.TemporaryDirectory(suffix=suffix, prefix=prefix)


def cleanup(dirname):
    shutil.rmtree(dirname)


def remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def dir_base(infile):
    basename = os.path.splitext(os.path.basename(infile))[0]
    dirname = os.path.dirname(infile)
    return dirname, basename


def dirbase(infile):
    return os.path.splitext(infile)[0]


def mkdir(indir):
    if not os.path.exists(indir):
        os.makedirs(indir)
    return indir


def timeit(f):
    def wrap(*args, **kwargs):
        start_time = time.time()
        ret = f(*args, **kwargs)
        end_time = time.time()
        if end_time - start_time < 60:
            tsc.logg.info(f"Time needed for {f.__name__}: {end_time - start_time:.3f}")
        else:
            tsc.logg.info(
                f"Time needed for {f.__name__}:"
                f" {time.strftime('%Hh:%Mm:%Ss', time.gmtime(end_time - start_time))}"
            )
        return ret

    return wrap


def get_file(key):
    components = key.split("/")
    return pkg_resources.resource_filename(components[0], "/".join(components[1:]))


def split_mut(mut):
    alt = mut.split(".")[-1]
    try:
        ref = mut.split(".")[-2]
    except:
        return None, None, None, None, None, None
    pos = mut.split(".")[-3]
    chrom = mut.split(".")[-4]
    gene = None
    ens = None
    try:
        gene = mut.split(".chr")[0].split("_")[1]
        ens = mut.split(".chr")[0].split("_")[0]
    except:
        pass
    return ens, gene, chrom, pos, ref, alt


def flips_in_cells(adata, df_in, df_out):
    c01 = ((df_in == 0) & (df_out == 1)).sum(axis=1)
    c10 = ((df_in == 1) & (df_out == 0)).sum(axis=1)
    c31 = ((df_in == 3) & (df_out == 1)).sum(axis=1)
    c30 = ((df_in == 3) & (df_out == 0)).sum(axis=1)
    c01 = c01 / c01.sum()
    c10 = c10 / c10.sum()
    c31 = c31 / c31.sum()
    c30 = c30 / c30.sum()
    c01 = c01.rename("f_0_1")
    c10 = c10.rename("f_1_0")
    c31 = c31.rename("f_3_1")
    c30 = c30.rename("f_3_0")
    df = pd.DataFrame([c01, c10, c31, c30]).T
    df["f_0_1_color"] = "#E41A1C"
    df["f_1_0_color"] = "#FF7F00"
    df["f_3_1_color"] = "#FDBF6F"
    df["f_3_0_color"] = "#FB9A99"
    adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True)


def get_database(species):
    data = []
    temp = {}
    with open(f"/home/frashidi/database/{species}/ensembl/{species}.rsem") as fin:
        i = -2
        for line in fin:
            i += 1
            line = line.strip()
            if i == -1:
                continue
            if i % 6 == 0:
                temp = {}
                temp["trans_ens"] = line.split("\t")[0]
                temp["trans_sym"] = line.split("\t")[1]
            if i % 6 == 1:
                temp["gene_ens"] = line.split("\t")[0]
                temp["gene_sym"] = line.split("\t")[1]
            if i % 6 == 2:
                temp["chrom"] = line
            if i % 6 == 3:
                temp["strand"] = line.split(" ")[0]
            if i % 6 == 4:
                temp["n_exons"] = int(line.split(" ")[0])
                temp["regions"] = line.split(" ")[1:]
                temp["start"] = int(line.split(" ")[1])
                temp["end"] = int(line.split(" ")[-1])
            if i % 6 == 5:
                details = line.split("; ")
                for det in details:
                    temp[det.split(" ")[0]] = det.split(" ")[1]
                data.append(temp)
    database = pd.DataFrame(data)
    return database


def translate_gene(alist, kind):
    mm10 = get_database("mm10")
    mm10 = mm10.set_index("gene_sym")["gene_ens"].drop_duplicates().to_dict()
    hg19 = get_database("hg19")
    hg19 = hg19.set_index("gene_sym")["gene_ens"].drop_duplicates().to_dict()
    if kind == "mm10":
        translator = mm10
    else:
        translator = hg19
    result = []
    for d in alist:
        if d in translator:
            result.append(f"{translator[d]}_{d}")
        else:
            tsc.logg.info(d)
    return result


def difference_between(dna, rna):
    dna_temp = (
        dna.var["chrom"].astype(str)
        + ":"
        + dna.var["position"].astype(str)
        + ":"
        + dna.var["reference"].astype(str)
        + ":"
        + dna.var["alteration"].astype(str)
    )
    rna_temp = (
        rna.var["chrom"].astype(str)
        + ":"
        + rna.var["position"].astype(str)
        + ":"
        + rna.var["reference"].astype(str)
        + ":"
        + rna.var["alteration"].astype(str)
    )
    a = np.intersect1d(dna_temp.values, rna_temp.values)
    b = np.setdiff1d(dna_temp.values, rna_temp.values)
    c = np.setdiff1d(rna_temp.values, dna_temp.values)
    tsc.logg.info(f"{len(b)}, inter={len(a)}, {len(c)}")

    b = dna.to_df(layer="mutant").T.loc[dna_temp[dna_temp.isin(b)].index]
    c = rna.to_df(layer="mutant").T.loc[rna_temp[rna_temp.isin(c)].index]

    return dna_temp[dna_temp.isin(a)].index, b, c


def with_timeout(timeout):
    def decorator(decorated):
        @functools.wraps(decorated)
        def inner(*args, **kwargs):
            pool = multiprocessing.pool.ThreadPool(1)
            async_result = pool.apply_async(decorated, args, kwargs)
            try:
                return async_result.get(timeout)
            except multiprocessing.TimeoutError:
                return None

        return inner

    return decorator


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
