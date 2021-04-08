import os

import pytest
from _helpers import *
from click.testing import CliRunner

import trisicell as tsc
from trisicell.commands.trisicell import cli


class TestCommands:
    def test_scistree(self):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "scistree",
                f"{tsc.ul.get_file('trisicell.datasets/test/test.SC')}",
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_scite(self):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "scite",
                f"{tsc.ul.get_file('trisicell.datasets/test/test.SC')}",
                "0.0000001",
                "0.1",
                "-r 3",
                "-l 1000",
            ],
        )
        assert result.exit_code == 0

    def test_phiscsb(self):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "phiscsb",
                f"{tsc.ul.get_file('trisicell.datasets/test/test.SC')}",
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_consensus(self):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "consensus",
                f"{tsc.ul.get_file('trisicell.datasets/test/consensus/biorxiv.fig3b.CFMatrix')}",
                f"{tsc.ul.get_file('trisicell.datasets/test/consensus/biorxiv.figs18a.CFMatrix')}",
                f"{tsc.ul.get_file('trisicell.datasets/test/consensus.CFMatrix')}",
            ],
        )
        assert result.exit_code == 0


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    def remove_test_dir():
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.scistree.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.scistree.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.scite.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.scite.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/consensus.CFMatrix"))

    request.addfinalizer(remove_test_dir)
