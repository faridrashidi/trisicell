import pytest
from _helpers import *

import trisicell as tsc


class TestCNA:
    @skip_mpi4py
    def test_infercna(self):
        assert 1 == 1
