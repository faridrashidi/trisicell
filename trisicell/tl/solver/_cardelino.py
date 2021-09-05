import trisicell as tsc


def cardelino(adata, mode, n_clones=None, tree=None):
    """Running Cardelino

    Computational integration of somatic clonal substructure and single-cell
    transcriptomes :cite:`Cardelino`.

    This method uses the prebuilt tree (either by bWES or scDNA) and
    then tries to cluster single-cells of scRNA based on the appearance
    of mutations detect in the tree.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of mutant and total.
    mode : :obj:`str`
        Cardelino mode, possible values are:

            - `fixed`: in this mode provide `tree`
            - `free`: in this mode provide `n_clones`
    n_clones : :obj:`int`, optional
        Number of clones, by default None
    tree : :class:`pandas.DataFrame`, optional
        Input prebuilt tree in the form of conflict-free matrix, by default None

    Returns
    -------
    :class:`pandas.DataFrame`
        The assignment of single-cells to the clones of the tree.
    """

    cardelino, cardelino_is_not_imported = tsc.ul.import_rpy2(
        "cardelino",
        "devtools::install_github('single-cell-genetics/cardelino'\n",
    )
    if cardelino_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    tsc.logg.info(f"running Cardelino")

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        A = ro.conversion.py2rpy(adata.to_df(layer="mutant").T)
        D = ro.conversion.py2rpy(adata.to_df(layer="total").T)
        if mode == "fixed":
            config = ro.conversion.py2rpy(tree.T)

    ro.globalenv["A"] = A
    ro.globalenv["D"] = D
    if mode == "fixed":
        ro.globalenv["config"] = config

    if mode == "free":
        assignments = cardelino.clone_id(A, D, n_clone=n_clones)
    elif mode == "fixed":
        assignments = cardelino.clone_id(A, D, Config=config)
    else:
        tsc.logg.error("Wrong mode")

    result = cardelino.assign_cells_to_clones(assignments.rx2["prob"])

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        result = ro.conversion.rpy2py(result)

    return result
