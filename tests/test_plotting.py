import pytest
from _helpers import *

import trisicell as tsc


class TestTrees:
    def __init__(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        self.tree = tsc.ul.to_tree(data)

    def test_clonal_tree(self):
        tsc.pl.clonal_tree(self.tree)

    @skip_mpi4py
    def test_dendro_tree(self):
        tsc.pl.dendro_tree(self.tree)
