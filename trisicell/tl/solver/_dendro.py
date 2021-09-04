import os

from IPython.display import Image, display

import trisicell as tsc


def dendro(adata):
    """Running DENDRO

    Genetic heterogeneity profiling and subclone detection by single-cell
    RNA sequencing :cite:`DENDRO`.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of `genotype`, `mutant` and `total`.

    Returns
    -------
    :obj:`None`
    """

    dendro, dendro_is_not_imported = tsc.ul.import_rpy2(
        "DENDRO", "devtools::install_github('zhouzilu/DENDRO')\n"
    )
    if dendro_is_not_imported:
        tsc.logg.error("Unable to import a package!")
    ggtree, ggtree_is_not_imported = tsc.ul.import_rpy2(
        "ggtree",
        "devtools::install_github(c('YuLab-SMU/ggtree','xiangpin/ggtreeExtra'"
        + ",'YuLab-SMU/aplot'))\ninstall.packages('cowplot')\n",
    )
    if ggtree_is_not_imported:
        tsc.logg.error("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.lib import grdevices
    from rpy2.robjects.packages import importr

    stats = importr("stats")

    tsc.logg.info(f"running DENDRO")

    Z = adata.to_df(layer="genotype")
    Z[(Z == 1) | (Z == 3)] = 1
    Z[Z == 2] = 0
    X = adata.to_df(layer="mutant")
    N = adata.to_df(layer="total")

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        Z = ro.conversion.py2rpy(Z.T)
        X = ro.conversion.py2rpy(X.T)
        N = ro.conversion.py2rpy(N.T)

    dist = dendro.DENDRO_dist(X, N, Z, show_progress=False)
    tree = stats.hclust(dist, method="ward.D")
    ro.globalenv["tree"] = tree

    cmd = f"""
    p <- ggtree(tree, layout='fan')
    p <- p + geom_tiplab()
    """

    with grdevices.render_to_bytesio(
        grdevices.png, width=800, height=800, res=100
    ) as image:
        p = ro.r(cmd)
        ro.r.show(p)

    return display(Image(image.getvalue(), embed=True, retina=True))
