import numpy as np
import pandas as pd
from IPython.display import SVG, Image, display

import trisicell as tsc


def infercna(expr, ref_cells, genome="hg19", show_heatmap=False):
    """Calling copy number profiles from scRNAseq.

    This function is a wraper to run
    `InferCNA <https://github.com/jlaffy/infercna>`_

    Parameters
    ----------
    expr : :class:`anndata.AnnData`
        The input expression matrix where

            - raw count are in `.X`
            - fpkm is in `.layer['fpkm']`
            - tpm is in `.layer['tpm']`
    ref_cells : :obj:`dict`
        A dictionary of list of normal cells.
    genome : :obj:`str`, optional
        {'hg38', 'hg19', 'mm10'}, by default "hg19"
    show_heatmap : :obj:`bool`, optional
        Show the copy number heatmap, by default False

    Returns
    -------
    :class:`pandas.DataFrame`
        A datafrme of log2 copy number ratios where normal cells
        are filtered out. Rows are cells and columns are list of
        variable genes.

    Examples
    --------
    >>> cna_log = tsc.tl.infercna(crc1, {'normal': ['NC_502']})
    """

    infercna, infercna_is_not_imported = tsc.ul.import_rpy2(
        "infercna",
        "devtools::install_github(c('jlaffy/scalop','jlaffy/infercna'))\n",
    )
    if infercna_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.lib import grdevices
    from rpy2.robjects.packages import importr

    base = importr("base")

    tsc.logg.info(f"running inferCNA")

    rc = {key: ro.StrVector(value) for key, value in ref_cells.items()}
    rc = ro.ListVector(rc)

    tpm = expr.to_df(layer="tpm").T
    # tpm.index = tpm.index.str.split("_").str[1]

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        data = ro.conversion.py2rpy(tpm)
    data = base.as_matrix(data)
    data.rownames = ro.StrVector(tpm.index.tolist())

    infercna.useGenome(genome)
    # ref = infercna.retrieveGenome()

    cna_mal = infercna.infercna(
        data,
        refCells=rc,
        window=100,
        n=5000,
        range=ro.FloatVector([-3, 3]),
        noise=0,
        center_method="median",
        isLog=False,
        verbose=False,
    )
    ro.globalenv["cna_mal"] = cna_mal
    ro.globalenv["rc"] = rc
    cna_mal = ro.r("cna_mal[, !colnames(cna_mal) %in% unlist(rc)]")

    cna_plot = infercna.cnaPlot(
        cna_mal,
        limits=ro.FloatVector([-0.5, 0.5]),
        cols=infercna.heatCols,
        ratio=0.5,
        x_name="Chromosome",
        y_name="Cell",
        legend_title="Inferred CNA\n[log2 ratio]",
        x_hide=ro.StrVector(["13", "18", "21", "Y"]),
        order_cells=False,
        subset_genes=ro.NULL,
        euclid_dist=False,
        angle=ro.NULL,
        x_angle=ro.NULL,
        y_angle=0,
        axis_rel=1,
        base_size=12,
        axis_title_size=12,
        axis_text_size=11,
        base_col="#073642",
        title=ro.NULL,
        subtitle=ro.NULL,
        caption=ro.NULL,
        text_size=12,
        y_hide=ro.NULL,
        tile_size=0.1,
        tile_col=ro.NULL,
        legend_position="right",
        legend_height=2,
        legend_width=0.6,
        legend_rel=0.9,
        legend_colour="black",
        legend_breaks=ro.NULL,
        legend_labels=ro.NULL,
        legend_justification="top",
        legend_title_position="bottom",
        legend_title_angle=ro.NULL,
        legend_title_rel=0.9,
    )

    # clones = infercna.findClones(
    #     cna_mal,
    #     prob=0.95,
    #     coverage=0.8,
    #     mode_size=10,
    #     clone_size=3,
    #     by="arm",
    #     bySampling=False,
    #     nsamp=2000,
    #     force_tries=False,
    #     verbose=False,
    # )
    # clones = {key: clones.rx2(key) for key in clones.names}

    # cor = infercna.refCorrect(cna=cna_mal, noise=ro.NULL, isLog=True)
    # print(cor)

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        cna_log = ro.conversion.rpy2py(cna_plot.rx2("data"))
    cna_log = pd.pivot(cna_log, index="Cell", columns="Gene", values="CNA")

    if show_heatmap:
        with ro.lib.grdevices.render_to_bytesio(
            grdevices.png, width=1024, height=896, res=150
        ) as image:
            ro.r.show(cna_plot.rx2("p"))
        display(Image(data=image.getvalue(), embed=True, retina=True))

    return cna_log
