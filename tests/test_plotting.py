import trisicell as tsc

from ._helpers import skip_graphviz, skip_rpy2


class TestTrees:
    @skip_graphviz
    def test_clonal_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.clonal_tree(tree)

    @skip_rpy2
    def test_dendro_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.dendro_tree(tree)
