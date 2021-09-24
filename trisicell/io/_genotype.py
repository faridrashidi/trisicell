import os

import anndata as ad
import pandas as pd

import trisicell as tsc


def read(filepath):
    """Read genotype matrix and read-count matrix.

    The genotype matrix must be in the in format of :class:`pandas.DataFrame`
    The read-count matrix must be in the format of :class:`anndata.AnnData`.

    Parameters
    ----------
    filepath : :obj:`str`
        The path to the file. The extension must be one of
        [`.SC`, `.CFMatrix`, `.h5ad`, `.h5ad.gz`]

    Returns
    -------
    :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        Depends on the format of the input file the output type is different.
    """

    ext = os.path.splitext(filepath)[-1]
    if ext in [".SC", ".CFMatrix", ".before_FP_FN_NA", ".tsv"]:
        sc = pd.read_table(filepath, index_col=0)
        if len(sc.columns) != len(set(sc.columns)):
            tsc.logg.error("Mutation ids must be unique!")
        return sc
    elif ext in [".h5ad", ".gz"]:
        return ad.read(filepath)
    else:
        tsc.logg.error("Extension is wrong!")


def write(obj, filepath):
    """Write genotype matrix or read-count matrix into a file.

    Parameters
    ----------
    obj : :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        The input object which is going to be written in a file.
    filepath : :obj:`str`
        The file path where the `obj` must be written in.
    """

    if isinstance(obj, pd.DataFrame):
        obj.index.name = "cellIDxmutID"
        obj.to_csv(filepath, sep="\t")
    elif isinstance(obj, ad.AnnData):
        obj.write(filepath + ".h5ad.gz", compression="gzip")
    else:
        tsc.logg.error("Object instance is wrong!")
