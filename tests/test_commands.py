import os

import pytest
from click.testing import CliRunner

import trisicell as tsc
from trisicell.commands.trisicell import cli


class TestSolvers:
    @pytest.mark.skip(reason="Needs log and cfmatrix to be removed.")
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
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.scistree.CFMatrix"))
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.scistree.log"))
        assert result.exit_code == 0

    @pytest.mark.skip(reason="Needs log and cfmatrix to be removed.")
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
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.scite.CFMatrix"))
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.scite.log"))
        assert result.exit_code == 0

    @pytest.mark.skip(reason="Needs log and cfmatrix to be removed.")
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
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.CFMatrix"))
        os.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.log"))
        assert result.exit_code == 0
