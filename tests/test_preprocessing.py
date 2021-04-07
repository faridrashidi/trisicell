import pytest
from _helpers import *

import trisicell as tsc


class TestPreProcessing:
    @skip_gurobi
    def test_bifilter_1(self):
        df_in = tsc.datasets.test()
        df_filtered = tsc.pp.bifilter(df_in, 0.5, 0.2)
        assert df_filtered.shape == (10, 4)
