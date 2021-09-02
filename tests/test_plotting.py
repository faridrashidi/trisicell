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

    @skip_rpy2
    def test_dendro_tree(self):
        file = tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
        data = tsc.io.read(file)
        tree = tsc.ul.to_tree(data)
        tsc.pl.dendro_tree(tree)
