import pandas as pd
from IPython.display import Image, display

import trisicell as tsc


def dendro(
    adata,
    cell_info=None,
    label_color="group_color",
    layout="fan",
    tiplab_size=2.5,
    width=800,
    height=800,
    dpi=100,
):
    """Running DENDRO

    Genetic heterogeneity profiling and subclone detection by single-cell
    RNA sequencing :cite:`DENDRO`.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of `genotype`, `mutant` and `total`.
    cell_info : :class:`pandas.DataFrame`, optional
        Information about cells such as color, by default None
    label_color : :obj:`str`, optional
        The column name in which colors of cells are stored
        in the dataframe provided as `cell_info`, by default "group_color"
    tiplab_size : :obj:`float`, optional
        Cell name size in the tree, by default 2.5
    layout : :obj:`str`, optional
        Layout of the tree, by default "fan"
        Values are:

            - `fan`
            - `dendrogram`
    width : :obj:`int`, optional
        Width of the figure, by default 800
    height : :obj:`int`, optional
        Height of the figure, by default 800
    dpi : :obj:`int`, optional
        The resolution, by default 100

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

        if cell_info is not None:
            cell_info_r = ro.conversion.py2rpy(cell_info.reset_index())
        else:
            cell_info_r = ro.conversion.py2rpy(
                pd.DataFrame.from_dict(
                    {
                        "x": adata.obs.index,
                        "group_color": adata.obs.shape[0] * ["#000000"],
                    }
                )
            )
        ro.globalenv["info"] = cell_info_r

    dist = dendro.DENDRO_dist(X, N, Z, show_progress=False)
    tree = stats.hclust(dist, method="ward.D")
    ro.globalenv["tree"] = tree

    cmd = f"""
    p <- ggtree(tree, layout='{layout}', size=0.4) %<+% info
    p <- p + geom_tiplab(
        aes(color={label_color}), align=TRUE, linesize=.1, angle=-90,
        hjust=-0.1, size={tiplab_size}, fontface='bold'
    )
    p <- p + scale_color_identity(guide='none')
    """

    with grdevices.render_to_bytesio(
        grdevices.png, width=width, height=height, res=dpi
    ) as image:
        p = ro.r(cmd)
        ro.r.show(p)

    return display(Image(image.getvalue(), embed=True, retina=True))
