import trisicell as tsc


class TestParitionFunction:
    def test_partition_function(self):
        df_in = tsc.datasets.test()
        probs = tsc.tl.partition_function(
            df_in,
            alpha=0.000001,
            beta=0.1,
            n_samples=100,
            n_batches=10,
            muts=["mut12"],
            cells=["cell6", "cell17"],
        )
        assert probs.mean(axis=1).round(4).values[0] >= 0
