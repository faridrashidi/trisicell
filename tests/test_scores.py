import numpy as np

import trisicell as tsc


class TestScores:
    def test_ad(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        grnd = tsc.io.read(f1)
        sol = tsc.io.read(f2)
        ad = tsc.tl.ad(grnd, sol)
        assert np.abs(ad - 0.9778) < 0.0001

    def test_dl(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        grnd = tsc.io.read(f1)
        sol = tsc.io.read(f2)
        dl = tsc.tl.dl(grnd, sol)
        assert np.abs(dl - 0.9880) < 0.0001

    def test_mltd(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        grnd = tsc.io.read(f1)
        sol = tsc.io.read(f2)
        mltd = tsc.tl.mltd(grnd, sol)
        assert np.abs(mltd["normalized_similarity"] - 0.7800) < 0.0001

    def test_tpted(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
        grnd = tsc.io.read(f1)
        sol = tsc.io.read(f2)
        tpted = tsc.tl.tpted(grnd, sol)
        assert tpted == 0.99
