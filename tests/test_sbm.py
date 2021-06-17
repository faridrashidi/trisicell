import trisicell as tsc

from ._helpers import skip_graph_tool


class TestSBM:
    @skip_graph_tool
    def test_sbm(self):
        data = tsc.datasets.test()
        out = tsc.tl.sbm(data)
        tree = tsc.ul.to_tree(out)
        assert len(tree.nodes) == 3
