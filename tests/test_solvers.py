import pytest
from _helpers import *

import trisicell as tsc


class TestSolvers:
    @skip_rpy2
    def test_simulate(self):
        df_in = tsc.datasets.simulate(
            n_cells=10, n_muts=10, n_clones=3, alpha=0.00001, beta=0.2, missing=0.1
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_in)
        assert is_cf == False

    def test_scistree(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.scistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == True

    def test_scite(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.scite(
            df_in, alpha=0.0000001, beta=0.1, n_restarts=3, n_iters=1000
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == True

    @skip_mpi4py
    def test_bnb(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.bnb(df_in, bounding="simulated")
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == True

    def test_phiscsb(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == True

    @skip_gurobi
    def test_phiscs_original(self):
        adata = tsc.datasets.acute_lymphocytic_leukemia2()
        adata.var["VAF"] = (
            2
            * adata.var["MutantCount"]
            / (adata.var["MutantCount"] + adata.var["ReferenceCount"])
        )
        df_out = tsc.tl.phiscs_original(
            adata.to_df(),
            alpha=0.001,
            beta=0.181749,
            delta=0.2,
            kmax=3,
            kel_weight=0,
            vaf_info=adata.var[["VAF"]],
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == True

    def test_booster(self):
        df_in = tsc.datasets.test()
        df_out = tsc.tl.booster(
            df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="SCITE",
            sample_on="muts",
            sample_size=10,
            n_samples=20,
            begin_sample=0,
            n_jobs=1,
            time_out=10000,
            save_inter=False,
            dir_inter=".",
            base_inter=None,
            disable_tqdm=False,
            weight=5,
            no_subsampling=False,
            no_dependencies=False,
        )
        is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
        assert is_cf == False
