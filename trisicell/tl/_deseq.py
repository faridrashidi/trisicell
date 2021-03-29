import trisicell as tsc


def deseq():
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
    except:
        raise SystemExit(
            "A problem was encountered importing `rpy2`. "
            "To run this `rpy2` and `R` need to be installed."
        )
    try:
        importr("Seurat")
    except:
        raise SystemExit(
            "A problem was encountered importing `Seurat` in R. "
            "To run this `Seurat` needs to be installed in R. "
            "Use the following lines to installed them.\n\n"
            "BiocManager::install('Seurat')\n"
        )

    cmd = f"""
    suppressPackageStartupMessages({{
        library(Seurat)
    }})
    """
    robjects.r(cmd)

    # https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
    # https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindMarkers
