import trisicell as tsc


class TestSolversTmp:
    def test_iscistree(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.iscistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_siclonefit(self):
        assert 1 == 1

    def test_infscite(self):
        assert 1 == 1
