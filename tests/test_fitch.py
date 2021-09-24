import trisicell as tsc


class TestFitchAlgorithm:
    def test_titch_algorithm(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.scite(
            df_in, alpha=0.00001, beta=0.1, n_restarts=3, n_iters=1000
        )
        tree = tsc.ul.to_tree(df_out)
        for n in tree.nodes:
            if tsc.ul.is_leaf(tree, n):
                tree.nodes[n]["profile"] = [1]
        tsc.tl.fitch(tree)
        assert True
