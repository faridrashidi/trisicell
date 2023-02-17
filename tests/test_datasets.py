import trisicell as tsc

from ._helpers import skip_rpy2


class TestDatasets:
    @skip_rpy2("oncoNEM")
    def test_simulate_1(self):
        df_in = tsc.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0.001, beta=0.4, missing=0.2
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_in)
        assert not is_cf

    @skip_rpy2("oncoNEM")
    def test_simulate_2(self):
        df_ground = tsc.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0, beta=0, missing=0
        )
        df_noisy = tsc.datasets.add_noise(df_ground, alpha=0.001, beta=0.4, missing=0.2)
        assert not tsc.ul.is_conflict_free_gusfield(df_noisy)

    def test_load_datasets(self):
        adata = tsc.datasets.example()
        assert adata.shape == (83, 452)
        adata = tsc.datasets.test()
        assert adata.shape == (20, 20)

        # adata = tsc.datasets.sublines_bwes()
        # assert adata.shape == (24, 6653)
        # mdata = tsc.datasets.sublines_bwts()
        # assert mdata.shape == (33, 55937)
        # mdata = tsc.datasets.sublines_scrnaseq()
        # assert mdata.shape == (175, 55851)
        # mdata = tsc.datasets.treated_actla4()
        # assert mdata.shape == (508, 58710)
        # mdata = tsc.datasets.treated_igg_ss2()
        # assert mdata.shape == (163, 56854)
        # mdata = tsc.datasets.treated_igg_sw()
        # assert mdata.shape == (163, 56854)

        adata = tsc.datasets.colorectal2()
        assert adata.shape == (78, 25)
        adata = tsc.datasets.high_grade_serous_ovarian_cancer_3celllines()
        assert adata.shape == (891, 14068)
