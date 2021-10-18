import numpy as np
import pandas as pd

import trisicell as tsc

from ._helpers import skip_graphviz


class TestUtils:
    def test_hclustering_1(self):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randint(0, 2, size=(10, 10)))
        clusters = tsc.ul.hclustering(df)
        assert clusters[6].value_counts()[5] == 3

    def test_hclustering_2(self):
        adata = tsc.datasets.example()
        clusters = tsc.ul.hclustering(adata.to_df(), metric="l1")
        assert len(clusters) == 81

    def test_hclustering_3(self):
        adata = tsc.datasets.example()
        clusters = tsc.ul.hclustering(adata.to_df(), metric="cosine")
        assert len(clusters) == 81

    def test_dist_dendro(self):
        tsc.settings.verbosity = 0
        adata = tsc.datasets.example()
        dist = tsc.ul.dist_dendro(adata)
        assert dist.sum() > 262456
        assert dist.sum() < 262457

    def test_tree_to_cfmatrix(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = tsc.ul.to_tree(df_out)
        df_out2 = tsc.ul.to_cfmatrix(tree)
        df_out = df_out.loc[df_out2.index, df_out2.columns].copy()
        pd.testing.assert_frame_equal(df_out, df_out2, check_dtype=False)

    def test_tree_to_mtree(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = tsc.ul.to_tree(df_out)
        tree = tsc.ul.to_mtree(tree)
        assert len(tree.nodes) == 13

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

    def test_general(self):
        # tsc.ul.executable("kDPFC", "SPhyR")
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        tsc.ul.dir_base(file)
        tsc.ul.dirbase(file)
        tsc.ul.get_param(
            "simNo_2-s_7-m_20-h_1-minVAF_0.1-ISAV_0-n_10-fp_0-fn_0.1-na_0-d_0-l_1000000"
            ".SC"
        )

        @tsc.ul.timeit
        def _test():
            return None

        _test()

        assert True
