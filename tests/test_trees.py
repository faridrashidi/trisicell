import trisicell as tsc

from ._helpers import skip_graphviz


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

    @skip_graphviz
    def test_cells_muts_rooted_at(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)

        tree = tsc.ul.to_tree(data)
        assert len(tree.nodes) == 10
        tsc.pl.clonal_tree(
            tree, show_id=True, muts_as_number=False, cells_as_number=False
        )

        tree = tsc.pp.collapse(tree)
        assert len(tree.nodes) == 6
        tsc.pl.clonal_tree(
            tree, show_id=True, muts_as_number=False, cells_as_number=False
        )

        res = tsc.ul.cells_rooted_at(tree, "[8]")
        assert res.shape[0] == 33

        res = tsc.ul.muts_rooted_at(tree, "[8]")
        assert res.shape[0] == 51
