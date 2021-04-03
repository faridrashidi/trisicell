import networkx as nx
import pytest

import trisicell as tsc


class TestConsensus:
    def test_consensus_1(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3b.CFMatrix")
        f2 = tsc.ul.get_file(
            "trisicell.datasets/test/consensus/biorxiv.figs18a.CFMatrix"
        )
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        cnt_tree1, cnt_tree2 = tsc.tl.consensus_run(sc1, sc2)
        ii = nx.is_isomorphic(cnt_tree1, cnt_tree2)  # result in biorxiv.figs18b
        assert ii == True

    def test_consensus_2(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/recomb.fig1a.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/recomb.fig1b.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        cnt_tree1, cnt_tree2 = tsc.tl.consensus_run(sc1, sc2)
        ii = nx.is_isomorphic(cnt_tree1, cnt_tree2)  # result in recomb.fig1c
        assert ii == True

    def test_consensus_3(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig4b.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig4c.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        cnt_tree1, cnt_tree2 = tsc.tl.consensus_run(sc1, sc2)
        ii = nx.is_isomorphic(cnt_tree1, cnt_tree2)  # result in biorxiv.fig4d
        assert ii == True

    def test_consensus_4(self):
        f1 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3b.CFMatrix")
        f2 = tsc.ul.get_file("trisicell.datasets/test/consensus/biorxiv.fig3c.CFMatrix")
        sc1 = tsc.io.read(f1)
        sc2 = tsc.io.read(f2)
        cnt_tree1, cnt_tree2 = tsc.tl.consensus_run(sc1, sc2)
        ii = nx.is_isomorphic(cnt_tree1, cnt_tree2)  # result in biorxiv.fig3f
        assert ii == True
