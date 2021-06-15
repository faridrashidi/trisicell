import trisicell as tsc


class TestConsensus:
    def test_consensus_1(self):
        # result in biorxiv.figs18b
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3b.CFMatrix")
        f2 = tsc.ul.get_file(
            "trisicell.datasets/test/consensus/biorxiv.figs18a.CFMatrix"
        )
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        final_tree = tsc.tl.consensus_tree(sc1, sc2)
        assert len(final_tree.nodes) == 19

    def test_consensus_2(self):
        # result in recomb.fig1c
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/recomb.fig1a.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/recomb.fig1b.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        final_tree = tsc.tl.consensus_tree(sc1, sc2)
        assert len(final_tree.nodes) == 21

    def test_consensus_3(self):
        # result in biorxiv.fig4d
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig4b.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig4c.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        final_tree = tsc.tl.consensus_tree(sc1, sc2)
        assert len(final_tree.nodes) == 34

    def test_consensus_4(self):
        # result in biorxiv.fig3f
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3b.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3c.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        final_tree = tsc.tl.consensus_tree(sc1, sc2)
        assert len(final_tree.nodes) == 19
