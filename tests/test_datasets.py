import trisicell as tsc


class TestDatasets:
    def test_load_datasets(self):
        adata = tsc.datasets.acute_lymphocytic_leukemia1()
        assert adata.shape == (111, 20)
        adata = tsc.datasets.acute_lymphocytic_leukemia2()
        assert adata.shape == (102, 16)
        # complete the rest

    def test_load_signatures(self):
        markers = tsc.datasets.get_markers()
        assert markers.shape == (51, 3)
        sig = tsc.datasets.get_signatures("mm10")
        assert "ENSMUSG00000002602.16_Axl" in sig["perez_lineage_markers"]
        sig = tsc.datasets.get_signatures("hg19")
        assert "ENSG00000167601.7_AXL" in sig["perez_lineage_markers"]
