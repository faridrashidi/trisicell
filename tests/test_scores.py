import numpy as np

import trisicell as tsc


class TestScores:
    def setup_method(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        self.grnd = tsc.io.read(f1)
        self.sol = tsc.io.read(f2)

    def test_ad(self):
        ad = tsc.tl.ad(self.grnd, self.sol)
        assert np.abs(ad - 0.9778) < 0.0001

    def test_dl(self):
        dl = tsc.tl.dl(self.grnd, self.sol)
        assert np.abs(dl - 0.9880) < 0.0001

    def test_mltd(self):
        mltd = tsc.tl.mltd(self.grnd, self.sol)
        assert np.abs(mltd["normalized_similarity"] - 0.7800) < 0.0001

    def test_tpted(self):
        tpted = tsc.tl.tpted(self.grnd, self.sol)
        assert tpted == 0.99
