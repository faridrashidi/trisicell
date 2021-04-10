import pytest
from _helpers import *

import trisicell as tsc


class TestTrees:
    def test_to_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        assert len(list(tree.nodes)) == 10
        assert len(list(tree.edges)) == 9

    def test_mtree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        mtree = tsc.ul.to_mtree(tree)
        assert len(mtree.nodes[8]["label"]) == 13
