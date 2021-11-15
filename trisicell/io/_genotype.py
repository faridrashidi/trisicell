import os

import anndata as ad
import ete3
import networkx as nx
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
        [`.SC`, `.CFMatrix`, `.h5ad`, `.h5ad.gz`, `.nwk`]

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
    elif ext in [".nwk"]:
        return _read_nwk(filepath)
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


def _read_nwk(filepath):
    tree = ete3.Tree(filepath, format=1)
    G = nx.DiGraph()
    node2id = {}
    i = 0
    for n in tree.traverse("postorder"):
        if n.name == "" or "Inner" in n.name:
            G.add_node(i, label="––")
        else:
            G.add_node(i, label=str(n.name))
        node2id[n] = i
        i += 1

    for p in tree.traverse("postorder"):
        pn = node2id[p]
        for c in p.children:
            cn = node2id[c]
            G.add_edge(pn, cn)

    i = 0
    for e, u, _ in G.edges.data("label"):
        G.edges[(e, u)]["label"] = f"m{i}"
        i += 1
    G.graph["normal_cells"] = []
    G.graph["splitter_mut"] = "\n"
    G.graph["splitter_cell"] = "\n"
    data = tsc.ul.to_cfmatrix(G)
    return data
