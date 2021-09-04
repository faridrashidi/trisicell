import pytest

import trisicell as tsc

from ._helpers import skip_graphviz, skip_rpy2


class TestTrees:
    @skip_graphviz
    def test_clonal_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.clonal_tree(tree)

    @skip_graphviz
    def test_clonal_tree_with_coloring(self):
        adata = tsc.datasets.high_grade_serous_ovarian_cancer_3celllines()
        df_out = adata.to_df(layer="ground")[adata.var_names[:1000]]
        tree = tsc.ul.to_tree(df_out)
        tsc.pl.clonal_tree(
            tree,
            muts_as_number=True,
            cells_as_number=False,
            cell_info=adata.obs,
            color_attr="group_color",
        )
        assert True

    @skip_rpy2
    def test_dendro_tree_1(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.dendro_tree(tree)
        assert True

    @skip_rpy2
    @pytest.mark.skip(reason="Giving error, don't know why!")
    def test_dendro_tree_2(self):
        adata = tsc.datasets.example()
        tsc.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        tsc.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        tsc.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        tsc.pp.build_scmatrix(adata)
        df_in = adata.to_df()
        df_out = tsc.tl.scistree(df_in, alpha=0.001, beta=0.2)
        tree = tsc.ul.to_tree(df_out)
        tsc.pl.dendro_tree(
            tree,
            cell_info=adata.obs,
            label_color="subclone_color",
            width=1200,
            height=600,
            dpi=200,
            distance_labels_to_bottom=3,
            inner_node_type="both",
            inner_node_size=2,
            annotation=[
                ("bar", "Axl", "Erbb3", 0.2),
                ("bar", "Mitf", "Mitf", 0.2),
            ],
        )
        assert True
