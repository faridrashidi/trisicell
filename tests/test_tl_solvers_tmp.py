import pytest

import trisicell as tsc

from ._helpers import skip_graph_tool, skip_rpy2


class TestSolversTmp:
    def test_titch_algorithm(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.scite(
            df_in, alpha=0.00001, beta=0.1, n_restarts=3, n_iters=1000
        )
        tree = tsc.ul.to_tree(df_out)
        for n in tree.nodes:
            if tsc.ul.is_leaf(tree, n):
                tree.nodes[n]["profile"] = [1]
        tsc.tl.fitch(tree)
        assert True

    def test_rscistree(self):
        adata = tsc.datasets.colorectal2(readcount=True)
        df_out = tsc.tl.rscistree(adata, mode="haploid")
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_iscistree(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.iscistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_siclonefit(self):
        assert True

    def test_infscite(self):
        assert True

    @skip_graph_tool
    def test_sbm(self):
        data = tsc.datasets.test()
        out = tsc.tl.sbm(data)
        tree = tsc.ul.to_tree(out)
        assert len(tree.nodes) == 3

    @skip_rpy2()
    @pytest.mark.skip(reason="Unable to import a package on GitHub!")
    def test_infercna(self):
        expr = tsc.datasets.example(is_expression=True)
        df_cna = tsc.tl.infercna(expr, ref_cells={"normal": ["C15_1"]}, genome="mm10")
        df_cna.loc["C15_1"] = 0
        expr.obsm["cna"] = df_cna.loc[expr.obs_names]
        tsc.pl.heatmap(expr, layer="cna")

    @skip_rpy2()
    @pytest.mark.skip(reason="Unable to import a package on GitHub!")
    def test_dendro(self):
        adata = tsc.datasets.example()
        tsc.tl.dendro(adata)
        assert True

    @skip_rpy2()
    @pytest.mark.skip(reason="Takes 6 minutes!")
    def test_cardelino(self):
        adata = tsc.datasets.example()
        tsc.tl.cardelino(adata, mode="free", n_clones=11)
        assert True
