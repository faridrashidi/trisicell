import numpy as np
import pandas as pd

import trisicell as tsc


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
