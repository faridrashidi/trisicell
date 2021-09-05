import trisicell as tsc

from ._helpers import skip_gurobi


class TestPreProcessing:
    @skip_gurobi
    def test_bifiltering(self):
        df_in = tsc.datasets.test()
        df_filtered = tsc.pp.bifiltering(df_in, 0.5, 0.2)
        assert df_filtered.shape == (10, 4)

    def test_preprocessing(self):
        tsc.settings.verbosity = 0
        adata = tsc.datasets.example()
        tsc.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        tsc.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        tsc.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        assert adata.shape == (83, 267)
