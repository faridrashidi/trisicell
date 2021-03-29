import trisicell as tsc


def infer_cna(adata):
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
    except:
        raise SystemExit(
            "A problem was encountered importing `rpy2`. "
            "To run this `rpy2` and `R` need to be installed."
        )
    try:
        importr("infercna")
    except:
        raise SystemExit(
            "A problem was encountered importing `infercna` in R. "
            "To run this `infercna` needs to be installed in R. "
            "Use the following lines to installed them.\n\n"
            "devtools::install_github(c('jlaffy/scalop','jlaffy/infercna'))\n"
        )

    tsc.logg.info(f"running InferCNA")
    tmpdir = tsc.ul.tmpdir(prefix="trisicell.", suffix=".cn", dirname=".")

    cmd = f"""
    suppressMessages(library(infercna))
    suppressMessages(library(data.table))
    suppressMessages(library(ggplot2))
    useGenome('hg19')
    ref <- retrieveGenome()

    tpm <- fread('~/crc1/expr.tpm.tsv') %>% as.data.frame()
    row.names(tpm) <- tpm$V1
    tpm$V1 <- NULL
    tpm <- t(as.matrix(tpm))
    tpm <- log2(tpm/10 + 1)
    row.names(tpm) <- sapply(row.names(tpm),
                             function(x) strsplit(x, '_')[[1]][2],
                             USE.NAMES=FALSE)

    refCells <- list('normalcolon' = c('crc1_NC_502'))

    cna = infercna(m=tpm, refCells=refCells, n=5000, noise=0.1,
                   isLog=TRUE, verbose=FALSE)
    cnaM = cna[, !colnames(cna) %in% unlist(refCells)]
    obj = cnaPlot(cna = cnaM)

    ggsave(obj$p, filename='crc1.png', device='png',
           width=7.5, height=6, units='in', limitsize=FALSE)
    """
    robjects.r(cmd)

    tsc.ul.cleanup(tmpdir)
