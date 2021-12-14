import numpy as np
import pytest

import trisicell as tsc


class TestScores:
    def setup_method(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        self.grnd = tsc.io.read(f1)
        self.sol = tsc.io.read(f2)

    def test_gs(self):
        gs = tsc.tl.gs(self.grnd, self.sol)
        assert np.abs(gs - 0.9895) < 0.0001

    def test_ad(self):
        ad = tsc.tl.ad(self.grnd, self.sol)
        assert np.abs(ad - 0.9778) < 0.0001

    def test_dl(self):
        dl = tsc.tl.dl(self.grnd, self.sol)
        assert np.abs(dl - 0.9880) < 0.0001

    def test_cc(self):
        tsc.tl.cc(self.grnd, self.sol)
        assert True

    def test_tpted(self):
        tpted = tsc.tl.tpted(self.grnd, self.sol)
        assert tpted == 0.99

    def test_caset(self):
        caset = tsc.tl.caset(self.grnd, self.sol)
        assert np.abs(caset - 0.7847) < 0.0001

    def test_disc(self):
        disc = tsc.tl.disc(self.grnd, self.sol)
        assert np.abs(disc - 0.7762) < 0.0001

    def test_mp3(self):
        mp3 = tsc.tl.mp3(self.grnd, self.sol)
        assert np.abs(mp3 - 0.6582) < 0.001

    def test_rf(self):
        rf = tsc.tl.rf(self.grnd, self.sol)
        assert np.abs(rf - 0.4864) < 0.0001

    @pytest.mark.skip(
        reason="Using MLTD in two tests is taking so long in test_scores!"
    )
    def test_mltd(self):
        mltd = tsc.tl.mltd(self.grnd, self.sol)
        assert np.abs(mltd["normalized_similarity"] - 0.7800) < 0.0001
