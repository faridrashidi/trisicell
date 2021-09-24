import trisicell as tsc

from ._helpers import skip_gurobi


class TestPreProcessing:
    @skip_gurobi
    def test_bifiltering(self):
        df_in = tsc.datasets.test()
        df_filtered = tsc.pp.bifiltering(df_in, 0.5, 0.2)
        assert df_filtered.shape == (10, 4)

    def test_binary(self):
        tsc.settings.verbosity = (
            0  # UnicodeEncodeError: 'ascii' codec can't encode character '\xd7'
        )
        df_in = tsc.datasets.test()
        tsc.pp.binarym_filter_private_mutations(df_in)
        tsc.pp.binarym_filter_clonal_mutations(df_in)
        tsc.pp.binarym_filter_nonsense_mutations(df_in)
        tsc.pp.binarym_statistics(df_in)
        assert df_in.shape == (20, 19)
        adata = tsc.datasets.example()
        tsc.pp.build_scmatrix(adata)
        df_in = tsc.pp.consensus_combine(adata.to_df())
        assert df_in.shape == (11, 452)

    def test_readcount(self):
        tsc.settings.verbosity = (
            0  # UnicodeEncodeError: 'ascii' codec can't encode character '\xd7'
        )
        adata = tsc.datasets.example()
        tsc.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        tsc.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        tsc.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        tsc.pp.statistics(adata)
        assert adata.shape == (83, 267)
