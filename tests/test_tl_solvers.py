import trisicell as tsc

from ._helpers import skip_gurobi, skip_mpi4py, skip_rpy2


class TestSolvers:
    def setup_method(self):
        self.df_in = tsc.datasets.test()

    def test_scite(self):
        df_out = tsc.tl.scite(
            self.df_in, alpha=0.0000001, beta=0.1, n_restarts=3, n_iters=1000
        )
        assert tsc.ul.is_conflict_free_gusfield(df_out)
        assert tsc.ul.is_conflict_free(df_out)

    @skip_mpi4py
    def test_bnb_simulated(self):
        df_out = tsc.tl.bnb(self.df_in, bounding="simulated")
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    @skip_mpi4py
    def test_bnb_real(self):
        df_out = tsc.tl.bnb(self.df_in, bounding="real", time_limit=20)
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    def test_scistree(self):
        df_out = tsc.tl.scistree(self.df_in, alpha=0.0000001, beta=0.1)
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    def test_phiscsb(self):
        df_out = tsc.tl.phiscsb(self.df_in, alpha=0.0000001, beta=0.1)
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    # def test_huntress_both(self):
    #     df_out = tsc.tl.huntress(self.df_in, alpha=0.0000001, beta=0.1)
    #     assert tsc.ul.is_conflict_free_gusfield(df_out)

    # def test_huntress_fn(self):
    #     df_out = tsc.tl.huntress(self.df_in, alpha=0, beta=0.1)
    #     assert tsc.ul.is_conflict_free_gusfield(df_out)

    @skip_rpy2("oncoNEM")
    def test_onconem(self):
        df_out = tsc.tl.onconem(self.df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    @skip_gurobi
    def test_phiscsi(self):
        df_out = tsc.tl.phiscsi(self.df_in, alpha=0.0000001, beta=0.1)
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    @skip_gurobi
    def test_phiscsi_bulk_1(self):
        adata = tsc.datasets.acute_lymphocytic_leukemia2()
        adata.var["VAF"] = (
            2
            * adata.var["MutantCount"]
            / (adata.var["MutantCount"] + adata.var["ReferenceCount"])
        )
        df_out = tsc.tl.phiscsi_bulk(
            adata.to_df(),
            alpha=0.001,
            beta=0.181749,
            delta=0.2,
            kmax=3,
            vaf_info=adata.var[["VAF"]],
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    @skip_gurobi
    def test_phiscsi_bulk_2(self):
        adata = tsc.datasets.colorectal2()
        df_in = adata.to_df()
        alpha = adata.uns["params_fig7a"]["alpha"]
        beta = adata.uns["params_fig7a"]["beta"]
        df_out = tsc.tl.phiscsi_bulk(df_in, alpha, beta, time_limit=120)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        flips_0_1, _, _, _ = tsc.ul.count_flips(df_in.values, df_out.values)
        assert is_cf
        assert flips_0_1 == 150

    @skip_gurobi
    def test_phiscsi_bulk_3(self):
        adata = tsc.datasets.colorectal2()
        df_in = adata.to_df()
        alpha = adata.uns["params_fig7b"]["alpha"]
        beta = adata.uns["params_fig7b"]["beta"]
        kmax = adata.uns["params_fig7b"]["kmax"]
        df_out = tsc.tl.phiscsi_bulk(df_in, alpha, beta, kmax, time_limit=120)
        assert df_out.columns[df_out.sum() == 0][0] == "ATP7B_chr13_52534322"

    @skip_gurobi
    def test_phiscs_readcount(self):
        adata = tsc.datasets.colorectal2(readcount=True)
        df_out = tsc.tl.phiscs_readcount(adata, alpha=0.01, beta=0.19)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_booster_phiscs(self):
        df_out = tsc.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="PhISCS",
            sample_on="muts",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            time_limit=120,
            dep_weight=5,
        )
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    def test_booster_scite(self):
        df_out = tsc.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="SCITE",
            sample_on="muts",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            n_iterations=10000,
            dep_weight=5,
        )
        assert tsc.ul.is_conflict_free_gusfield(df_out)

    def test_booster_scistree_on_cells(self):
        df_out = tsc.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="ScisTree",
            sample_on="cells",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            n_iterations=10000,
            dep_weight=5,
        )
        assert tsc.ul.is_conflict_free_gusfield(df_out)
