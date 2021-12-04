import pytest
from click.testing import CliRunner

import trisicell as tsc
from trisicell.commands.trisicell import cli

from ._helpers import skip_graphviz


class TestCommands:
    def setup_method(self):
        self.runner = CliRunner()

    def test_scistree(self):
        result = self.runner.invoke(
            cli,
            [
                "scistree",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_huntress_both(self):
        result = self.runner.invoke(
            cli,
            [
                "huntress",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_scite(self):
        result = self.runner.invoke(
            cli,
            [
                "scite",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
                "-r 3",
                "-l 1000",
            ],
        )
        assert result.exit_code == 0

    def test_scite_experiment(self):
        result = self.runner.invoke(
            cli,
            [
                "scite",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
                "-e",
                "-t 1",
                "-s 1",
            ],
        )
        assert result.exit_code == 0

    def test_phiscsb(self):
        result = self.runner.invoke(
            cli,
            [
                "phiscsb",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_consensus(self):
        path = "trisicell.datasets/test/consensus"
        result = self.runner.invoke(
            cli,
            [
                "consensus",
                tsc.ul.get_file(path + "/biorxiv.fig3b.CFMatrix"),
                tsc.ul.get_file(path + "/biorxiv.figs18a.CFMatrix"),
                tsc.ul.get_file("trisicell.datasets/test/consensus.CFMatrix"),
            ],
        )
        assert result.exit_code == 0

    def test_cf2newick(self):
        result = self.runner.invoke(
            cli,
            [
                "cf2newick",
                tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.CFMatrix"),
            ],
        )
        assert result.exit_code == 0

    def test_mcalling(self):
        result = self.runner.invoke(
            cli,
            [
                "mcalling",
                tsc.ul.get_file("trisicell.datasets/test/mcalling_config.yml"),
                "--test",
            ],
        )
        assert result.exit_code == 0

    def test_search(self):
        result = self.runner.invoke(
            cli,
            ["search", tsc.ul.get_file("trisicell.datasets/test/test.tsv"), "-p 2"],
        )
        assert result.exit_code == 0

    def test_score(self):
        result = self.runner.invoke(
            cli,
            [
                "score",
                tsc.ul.get_file(
                    "trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix"
                ),
                tsc.ul.get_file(
                    "trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix"
                ),
            ],
        )
        assert result.exit_code == 0

    @skip_graphviz
    def test_cf2tree(self):
        result = self.runner.invoke(
            cli,
            [
                "cf2tree",
                tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.CFMatrix"),
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="Error Don't know!")
    def test_partf(self):
        result = self.runner.invoke(
            cli,
            [
                "partf",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0001",
                "0.1",
                "--n_threads 2",
                "--n_samples 100",
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="PyTest issue with multithreading!")
    def test_booster(self):
        result = self.runner.invoke(
            cli,
            [
                "booster",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "0.0000001",
                "0.1",
                "--solver scite",
                "--n_samples 100",
                "--sample_size 15",
                "--n_jobs 2",
                "--n_iterations 10000",
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="pyBnB issue with CLI!")
    def test_bnb(self):
        result = self.runner.invoke(
            cli,
            [
                "bnb",
                tsc.ul.get_file("trisicell.datasets/test/test.tsv"),
                "-b simulated",
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
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.huntress.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.huntress.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.booster.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.booster.log"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/consensus.CFMatrix"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.info2"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.newick"))
        tsc.ul.remove(tsc.ul.get_file("trisicell.datasets/test/test.phiscsb.png"))
        tsc.ul.cleanup(tsc.ul.get_file("trisicell.datasets/test/_map"))
        tsc.ul.cleanup(tsc.ul.get_file("trisicell.datasets/test/_tmp"))
        tsc.ul.cleanup(tsc.ul.get_file("trisicell.datasets/test/test"))
        # tsc.ul.cleanup(tsc.ul.get_file("trisicell.datasets/test/test.partf.samples"))

    request.addfinalizer(remove_test_dir)
