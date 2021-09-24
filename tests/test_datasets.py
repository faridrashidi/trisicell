import trisicell as tsc

from ._helpers import skip_rpy2


class TestDatasets:
    @skip_rpy2
    def test_simulate_1(self):
        df_in = tsc.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0.001, beta=0.4, missing=0.2
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_in)
        assert not is_cf

    @skip_rpy2
    def test_simulate_2(self):
        df_ground = tsc.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0, beta=0, missing=0
        )
        df_noisy = tsc.datasets.add_noise(df_ground, alpha=0.001, beta=0.4, missing=0.2)
        assert not tsc.ul.is_conflict_free_gusfield(df_noisy)

    def test_load_datasets(self):
        # TODO: complete
        adata = tsc.datasets.acute_lymphocytic_leukemia1()
        assert adata.shape == (111, 20)
        adata = tsc.datasets.acute_lymphocytic_leukemia2()
        assert adata.shape == (102, 16)
        adata = tsc.datasets.acute_lymphocytic_leukemia3()
        adata = tsc.datasets.acute_lymphocytic_leukemia4()
        adata = tsc.datasets.acute_lymphocytic_leukemia5()
        adata = tsc.datasets.acute_lymphocytic_leukemia6()
        # adata = tsc.datasets.colorectal1()
        adata = tsc.datasets.colorectal2()
        # adata = tsc.datasets.colorectal3()
        adata = tsc.datasets.erbc()
        adata = tsc.datasets.example()
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer1()
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer2()
        # adata = tsc.datasets.high_grade_serous_ovarian_cancer3()
        adata = tsc.datasets.high_grade_serous_ovarian_cancer_3celllines()
        adata = tsc.datasets.melanoma20()
        adata = tsc.datasets.muscle_invasive_bladder()
        adata = tsc.datasets.myeloproliferative_neoplasms18()
        adata = tsc.datasets.myeloproliferative_neoplasms78()
        adata = tsc.datasets.myeloproliferative_neoplasms712()
        adata = tsc.datasets.oligodendroglioma_idh_mutated_tumor()
        adata = tsc.datasets.renal_cell_carcinoma()
        adata = tsc.datasets.test()
        adata = tsc.datasets.tnbc()

    def test_load_signatures(self):
        markers = tsc.datasets.get_markers()
        assert markers.shape == (51, 3)
        sig = tsc.datasets.get_signatures("mm10")
        assert "ENSMUSG00000002602.16_Axl" in sig["perez_lineage_markers"]
        sig = tsc.datasets.get_signatures("hg19")
        assert "ENSG00000167601.7_AXL" in sig["perez_lineage_markers"]
