import pytest
from _helpers import *

import trisicell as tsc


class TestTrees:
    @skip_graphviz
    def test_clonal_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.clonal_tree(tree)

    @skip_mpi4py
    def test_dendro_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.dendro_tree(tree)
