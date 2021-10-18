import trisicell as tsc


class TestIO:
    def test_read_newick(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/T00.nwk")
        f2 = tsc.ul.get_file("trisicell.datasets/test/T06.nwk")
        df1 = tsc.io.read(f1)
        df2 = tsc.io.read(f2)
        assert df1.shape == (185, 368)
        assert df2.shape == (185, 368)
