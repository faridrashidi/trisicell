import pytest

import trisicell as tsc


class TestSolvers:
    def test_scistree(self):
        x = tsc.ul.add(1, 2)
        assert x == 3
