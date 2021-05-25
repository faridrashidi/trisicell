from _helpers import *

import trisicell as tsc


class TestSBM:
    @skip_graph_tool
    def test_sbm(self):
        data = tsc.datasets.test()
        out = tsc.tl.sbm(data)
        tree = tsc.ul.to_tree(out)
        assert len(tree.nodes) == 3
