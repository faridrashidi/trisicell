import os

import trisicell as tsc


def dendro(adata, kind):
    """Running DENDRO

    https://doi.org/10.1186/s13059-019-1922-x

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of mutant and total.
    kind : :obj:`str`
        Which dataset is this, either `brna` or `scrna`
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
        importr("DENDRO")
    except:
        raise SystemExit(
            "A problem was encountered importing `DENDRO` in R. "
            "To run this `DENDRO` needs to be installed in R. "
            "Use the following lines to installed them.\n\n"
            "devtools::install_github('zhouzilu/DENDRO')\n"
        )

    tsc.logg.info(f"running DENDRO")
    tmpdir = tsc.ul.tmpdir(suffix=".dendro", dirname=".")

    Z = adata.to_df(layer="genotype")
    Z[(Z == 1) | (Z == 3)] = 1
    Z[Z == 2] = 0
    X = adata.to_df(layer="mutant")
    N = adata.to_df(layer="total")
    Z.to_csv(os.path.join(tmpdir, "Z.csv"))
    X.to_csv(os.path.join(tmpdir, "X.csv"))
    N.to_csv(os.path.join(tmpdir, "N.csv"))

    if kind == "brna":
        width = 3.5
        height = 2
        body1 = f"""
        file <- './trisicell_notebooks/data/mouse.melanoma/24.bulk.rna.lcbg.csv'
        colors <- c('#D92347','#FF7F00','#812C7C','#009949','#0067A5')
        """

        body2 = f"""
        p <- ggtree(tree, layout='dendrogram', size=0.4) %<+% info1
        p <- p + geom_tiplab(aes(color=group), align=TRUE,
            linesize=.1, angle=-90, hjust=-0.1, size=2.5, fontface='bold')
        p <- p + scale_color_manual(values=colors, guide='none')
        p <- p + xlim_tree(20)
        """
    elif kind == "scrna":
        width = 8
        height = 3
        body1 = f"""
        file <- './trisicell_notebooks/data/mouse.melanoma/24.sc.rna.lcbg.csv'
        """
        body2 = f"""
        p <- ggtree(tree, layout='dendrogram', size=0.4) %<+% info1
        p <- p + geom_tiplab(aes(color=subclone_color), align=TRUE,
            linesize=.1, angle=-90, hjust=-0.1, size=1.5, fontface='bold')
        p <- p + scale_color_identity(guide='none')
        p <- p + xlim_tree(5)
        """
    else:
        width = 5
        height = 2
        body1 = f"""
        """
        body2 = f"""
        """

    cmd = f"""
    suppressMessages(library(DENDRO))
    suppressMessages(library(data.table))

    N <- t(as.matrix(fread('{tmpdir}/N.csv'), rownames=1))
    X <- t(as.matrix(fread('{tmpdir}/X.csv'), rownames=1))
    Z <- t(as.matrix(fread('{tmpdir}/Z.csv'), rownames=1))

    dist = DENDRO.dist(X, N, Z, show.progress=FALSE)

    # pdf('DENDRO.plot.pdf', width=8, height=8)
    # cluster = DENDRO.cluster(dist, label=colnames(Z), type='fan')
    # dev.off()

    suppressMessages(library(ggtree))
    suppressMessages(library(ggtreeExtra))
    suppressMessages(library(aplot))
    suppressMessages(library(ggplot2))
    suppressMessages(library(cowplot))

    {body1}

    info1 <- fread(file)
    info1 <- info1 %>% as.data.frame()
    tree <- hclust(dist, method='ward.D')

    # p <- ggtree(tree, layout='fan', size=0.5) %<+% info1 +
    #     geom_tippoint(mapping=aes(color=MPS), size=5) +
    #     scale_color_gradient2(low='#E41A1C', midpoint=0, mid='white', high='#377EB8')

    {body2}

    ggsave(p, filename='DENDRO.final.pdf', device='pdf',
        width={width}, height={height}, units='in', limitsize=FALSE)
    """
    robjects.r(cmd)

    tsc.ul.cleanup(tmpdir)
