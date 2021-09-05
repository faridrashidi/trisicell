import trisicell as tsc

from ._helpers import skip_graph_tool, skip_rpy2


class TestSolversTmp:
    def test_iscistree(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.iscistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_rscistree(self):
        # adata = tsc.datasets.colorectal2(readcount=True)
        # df_out = tsc.tl.rscistree(adata, mode="haploid")
        # is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        # assert is_cf
        assert True

    def test_siclonefit(self):
        assert True

    def test_infscite(self):
        assert True

    @skip_rpy2
    def test_dendro(self):
        adata = tsc.datasets.example()
        tsc.tl.dendro(adata)
        assert True

    @skip_graph_tool
    def test_sbm(self):
        data = tsc.datasets.test()
        out = tsc.tl.sbm(data)
        tree = tsc.ul.to_tree(out)
        assert len(tree.nodes) == 3
