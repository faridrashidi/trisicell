import trisicell as tsc

from ._helpers import skip_gurobi


class TestPreProcessing:
    @skip_gurobi
    def test_bifiltering(self):
        df_in = tsc.datasets.test()
        df_filtered = tsc.pp.bifiltering(df_in, 0.5, 0.2)
        assert df_filtered.shape == (10, 4)
