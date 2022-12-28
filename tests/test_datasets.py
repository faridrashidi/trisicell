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

        adata = tsc.datasets.acute_lymphocytic_leukemia1()
        assert adata.shape == (111, 20)
        adata = tsc.datasets.acute_lymphocytic_leukemia2()
        assert adata.shape == (102, 16)
        adata = tsc.datasets.acute_lymphocytic_leukemia3()
        assert adata.shape == (150, 49)
        adata = tsc.datasets.acute_lymphocytic_leukemia4()
        assert adata.shape == (143, 78)
        adata = tsc.datasets.acute_lymphocytic_leukemia5()
        assert adata.shape == (96, 105)
        adata = tsc.datasets.acute_lymphocytic_leukemia6()
        assert adata.shape == (146, 10)
        adata = tsc.datasets.colorectal1()
        assert adata.shape == (178, 16)
        adata = tsc.datasets.colorectal2()
        assert adata.shape == (78, 25)
        # adata = tsc.datasets.colorectal3()
        adata = tsc.datasets.erbc()
        assert adata.shape == (47, 40)
        adata = tsc.datasets.high_grade_serous_ovarian_cancer_3celllines()
        assert adata.shape == (891, 14068)
        adata = tsc.datasets.melanoma20()
        assert adata.shape == (20, 2367)
        adata = tsc.datasets.muscle_invasive_bladder()
        assert adata.shape == (44, 443)
        adata = tsc.datasets.myeloproliferative_neoplasms18()
        assert adata.shape == (58, 18)
        adata = tsc.datasets.myeloproliferative_neoplasms78()
        assert adata.shape == (58, 78)
        adata = tsc.datasets.myeloproliferative_neoplasms712()
        assert adata.shape == (58, 712)
        adata = tsc.datasets.oligodendroglioma_idh_mutated_tumor()
        assert adata.shape == (579, 77)
        adata = tsc.datasets.renal_cell_carcinoma()
        assert adata.shape == (17, 35)
        adata = tsc.datasets.tnbc()
        assert adata.shape == (16, 20)
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer1()
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer2()
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer3()
        # adata = tsc.datasets.acute_lymphocytic_leukemia_many()
