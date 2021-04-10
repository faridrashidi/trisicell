import numpy as np
import pandas as pd
from IPython.display import SVG, Image, display

import trisicell as tsc


def infercna(expr, ref_cells):
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

    tpm = expr.to_df(layer="tpm")
    tpm = np.log(tpm / 10 + 1).T
    tpm.index = tpm.index.str.split("_").str[1]

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        data = ro.conversion.py2rpy(tpm)
    data = base.as_matrix(data)
    data.rownames = ro.StrVector(tpm.index.tolist())

    infercna.useGenome("hg19")
    ref = infercna.retrieveGenome()
    cna = infercna.infercna(
        m=data, refCells=rc, n=5000, noise=0.1, isLog=True, verbose=False
    )
    print(type(cna))
    # cna_m = cna[, !colnames(cna) %in% base.unlist(refCells)]
    obj = infercna.cnaPlot(cna=cna)

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        data = ro.conversion.rpy2py(obj.rx2("data"))
    data = pd.pivot(data, index="Cell", columns="Gene", values="CNA")

    with ro.lib.grdevices.render_to_bytesio(
        grdevices.png, width=1024, height=896, res=150
    ) as image:
        ro.r.show(obj.rx2("p"))

    display(Image(data=image.getvalue(), embed=True, retina=True))

    return data
