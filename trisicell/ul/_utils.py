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
    icf = is_conflict_free_gusfield(df_out)
    tsc.logg.info("output -- CF: ", end="")
    if icf:
        tsc.logg.info(icf, color="green")
    else:
        tsc.logg.info(icf, color="red")
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


def calc_nll_matrix(df_in, df_out, alpha, beta):
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
    return -objective


def stat(df_in, df_out, alpha, beta, running_time):
    log_input(df_in)
    log_output(df_out, running_time)
    log_flip(df_in, df_out)
    nll = calc_nll_matrix(df_in, df_out, alpha, beta)
    tsc.logg.info(f"score -- NLL: {nll}")


def get_param(filename):
    def _get_param_helper(param):
        try:
            value = basename.split(f"{param}_")[1]
            if "-" in value:
                value = value.split("-")[0]
            else:
                value = value.split(".")[0]
            return float(value) if "." in value else int(value)
        except IndexError:
            return None

    data = {}
    _, basename = dir_base(filename)
    for param in [
        "simNo",
        "s",
        "m",
        "h",
        "minVAF",
        "ISAV",
        "n",
        "fp",
        "fn",
        "na",
        "d",
        "l",
    ]:
        value = _get_param_helper(param)
        if value is not None:
            data[param] = value
    # data["s"] = int(basename.split("s_")[1].split("-")[0])
    # data["m"] = int(basename.split("-")[2].split("_")[1])
    # data["h"] = int(basename.split("-")[3].split("_")[1])
    # data["minVAF"] = float(basename.split("-")[4].split("_")[1])
    # data["ISAV"] = int(basename.split("-")[5].split("_")[1])
    # data["n"] = int(basename.split("-")[6].split("_")[1])
    # data["fp"] = float(basename.split("-")[7].split("_")[1])
    # data["fn"] = float(basename.split("-")[8].split("_")[1])
    # data["na"] = float(basename.split("-")[9].split("_")[1])
    # data["d"] = float(basename.split("-")[10].split("_")[1])
    # last = basename.split("-")[11]
    # if "." in last:
    #     data["l"] = int(last.split(".")[0].split("_")[1])
    # else:
    #     data["l"] = int(last.split("_")[1])
    return data


def count_flips(I_mtr, O_mtr, na_value=3):
    flips_0_1 = 0
    flips_1_0 = 0
    flips_na_0 = 0
    flips_na_1 = 0
    n, m = I_mtr.shape
    for i in range(n):
        for j in range(m):
            if I_mtr[i, j] == 0 and O_mtr[i, j] == 1:
                flips_0_1 += 1
            elif I_mtr[i, j] == 1 and O_mtr[i, j] == 0:
                flips_1_0 += 1
            elif I_mtr[i, j] == na_value and O_mtr[i, j] == 0:
                flips_na_0 += 1
            elif I_mtr[i, j] == na_value and O_mtr[i, j] == 1:
                flips_na_1 += 1
    return flips_0_1, flips_1_0, flips_na_0, flips_na_1


def infer_rates(I_mtr, O_mtr, na_value=3):
    flips_0_1, flips_1_0, flips_na_0, flips_na_1 = count_flips(I_mtr, O_mtr, na_value)
    fn_rate = flips_0_1 / ((O_mtr == 1) & (I_mtr != na_value)).sum()
    fp_rate = flips_1_0 / ((O_mtr == 0) & (I_mtr != na_value)).sum()
    na_rate = (flips_na_1 + flips_na_0) / I_mtr.size
    return fn_rate, fp_rate, na_rate


def is_conflict_free(df_in):
    D = df_in.astype(int).values
    if not np.array_equal(np.unique(D), [0, 1]):
        return False
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


def is_conflict_free_gusfield(df_in):
    """Check conflict-free criteria via Gusfield algorithm.

    This is an implementation of algorithm 1.1 in :cite:`Gusfield_1991`.

    The order of this algorithm is O(nm)
    where n is the number of cells and m is the number of mutations.

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.

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

    I_mtr = df_in.astype(int).values
    if not np.array_equal(np.unique(I_mtr), [0, 1]):
        return False

    def _sort_bin(a):
        b = np.transpose(a)
        b_view = np.ascontiguousarray(b).view(
            np.dtype((np.void, b.dtype.itemsize * b.shape[1]))
        )
        idx = np.argsort(b_view.ravel())[::-1]
        c = b[idx]
        return np.transpose(c), idx

    Ip = I_mtr.copy()
    O_mtr, _ = _sort_bin(Ip)
    Lij = np.zeros(O_mtr.shape, dtype=int)
    for i in range(O_mtr.shape[0]):
        maxK = 0
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                Lij[i, j] = maxK
                maxK = j + 1
    Lj = np.amax(Lij, axis=0)
    for i in range(O_mtr.shape[0]):
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                if Lij[i, j] != Lj[j]:
                    return False
    return True


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


def executable(binary, appname):
    executable = shutil.which(binary)
    if executable is None:
        if not os.path.exists(f"{tsc.settings.tools}/{binary}"):
            tsc.logg.error(
                f"Cannot find the binary file of {appname} with `{binary}` name!"
            )
        else:
            executable = f"{tsc.settings.tools}/{binary}"
    return executable


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


def split_mut(mut):
    try:
        ref = mut.split(".")[-2]
        pos = mut.split(".")[-3]
        chrom = mut.split(".")[-4]
        gene = mut.split(".chr")[0].split("_")[1]
        ens = mut.split(".chr")[0].split("_")[0]
        alt = mut.split(".")[-1]
        return ens, gene, chrom, pos, ref, alt
    except Exception:
        return None, None, None, None, None, None
