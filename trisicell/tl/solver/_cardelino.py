import trisicell as tsc


def cardelino(adata, mode, kind, tree=None):
    """running Cardelino

    This method uses the prebuilt tree (either by bWES or scDNA) and
    then tries to cluster single-cells of scRNA based on the appearance
    of mutations detect in the tree.
    https://doi.org/10.1038/s41592-020-0766-3

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of mutant and total.
    mode : :obj:`str`
        Cardelino mode, either `fixed` or `free`
    kind : :obj:`str`
        Which dataset is this, either `brna` or `scrna`
    tree : :class:`pandas.DataFrame`
        Input prebuilt tree in the form of conflict-free matrix.

    Returns
    -------
    :class:`pandas.DataFrame`
        The assignment of single-cells to the clones of the tree.
    """

    try:
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
    except:
        raise SystemExit(
            "A problem was encountered importing `rpy2`. "
            "To run this `rpy2` and `R` need to be installed."
        )
    try:
        importr("cardelino")
    except:
        raise SystemExit(
            "A problem was encountered importing `cardelino` in R. "
            "To run this `cardelino` needs to be installed in R. "
            "Use the following lines to installed them.\n\n"
            "devtools::install_github('single-cell-genetics/cardelino',build_vignettes = FALSE)\n"
        )

    tsc.logg.info(f"running Cardelino")

    if kind == "brna":
        n_clone = 12
    elif kind == "scrna":
        n_clone = 25
    else:
        n_clone = 0

    A = adata.to_df(layer="mutant").T
    D = adata.to_df(layer="total").T
    if mode == "fixed":
        Config = tree.T

    robjects.pandas2ri.activate()
    A_r = robjects.conversion.py2rpy(A)
    D_r = robjects.conversion.py2rpy(D)
    if mode == "fixed":
        Config_r = robjects.conversion.py2rpy(Config)

    robjects.globalenv["A"] = A_r
    robjects.globalenv["D"] = D_r
    if mode == "fixed":
        robjects.globalenv["Config"] = Config_r

    cmd = f"""
    suppressPackageStartupMessages({{
        library(cardelino)
    }})
    A = as.matrix(A)
    D = as.matrix(D)
    """

    if mode == "free":
        cmd += f"""
        assignments <- clone_id(A, D, n_clone={n_clone})
        """
    elif mode == "fixed":
        cmd += f"""
        Config = as.matrix(Config)
        assignments <- clone_id(A, D, Config=Config)
        """
    else:
        raise ValueError("Wrong mode")

    cmd += f"""
    result <- assign_cells_to_clones(assignments$prob)
    result
    """

    # cmd = f"""
    # suppressPackageStartupMessages({{
    #     library(cardelino)
    # }})
    # vcf <- read_vcf(system.file('extdata', 'cell_example.mpileup.vcf.gz', package = 'cardelino'))
    # input_data <- get_snp_matrices(vcf)
    # canopy <- readRDS(system.file('extdata', 'canopy_results.example.rds', package = 'cardelino'))
    # Config <- canopy$tree$Z
    # rownames(Config) <- gsub('chr', '', gsub(':', '_', gsub('_.*', '', rownames(Config))))
    # assignments <- clone_id(input_data$A, input_data$D, Config) #mutant, total, tree
    # result <- assign_cells_to_clones(assignments$prob)
    # """

    result = robjects.r(cmd)

    return result
