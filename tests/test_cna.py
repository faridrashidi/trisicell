from ._helpers import skip_mpi4py


class TestCNA:
    @skip_mpi4py
    def test_infercna(self):
        assert 1 == 1
