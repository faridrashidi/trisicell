import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

import trisicell as tsc


def count_mut_in_a_gene(adata, cut_off=0):
    x = pd.Series([".".join(x.split(".")[:-4]) for x in adata.var_names]).value_counts()
    return x[x > cut_off]


def remove_indels(adata):
    snps = []
    for x in adata.var_names:
        ens, gene, chrom, pos, ref, alt = tsc.ul.split_mut(x)
        if ("-" not in ref) and ("-" not in alt):
            snps.append(x)
    keep_mut_by_list(adata, snps)


def most_frequent_muts(adata):
    G = adata.to_df(layer="genotype")
    V = adata.to_df(layer="mutant")
    T = adata.to_df(layer="total")
    s1 = ((G == 1) | (G == 3)).sum(axis=0).rename("n_cells")
    s2 = V.sum(axis=0).rename("overall_sum_mutant_coverage")
    s3 = T.sum(axis=0).rename("overall_sum_total_coverage")
    return pd.concat([s1, s2, s3], axis=1).sort_values(
        by=["n_cells", "overall_sum_mutant_coverage", "overall_sum_total_coverage"],
        ascending=False,
    )


def muts_in_a_subset_of_cells(adata, alist):
    G = adata[alist].layers["genotype"]
    return adata.var_names[((G == 1) | (G == 3)).sum(axis=0) > 0]


def n_muts_per_cell_with_vaf(adata):
    G = adata.layers["genotype"]
    n_muts = ((G == 1) | (G == 3)).sum(axis=1)
    V = adata.layers["mutant"]
    T = adata.layers["total"]

    vaf = np.zeros(G.shape)
    vaf = np.divide(V, T, where=T != 0)
    vaf[(G == 1) | (G == 3)] = 0
    vaf = np.nanmean(np.where(vaf != 0, vaf, np.nan), axis=1)

    n_mutant = V
    n_mutant[(G == 1) | (G == 3)] = 0
    n_mutant = np.nanmean(np.where(n_mutant != 0, n_mutant, np.nan), axis=1)

    n_total = T
    n_total[(G == 1) | (G == 3)] = 0
    n_total = np.nanmean(np.where(n_total != 0, n_total, np.nan), axis=1)
    return pd.DataFrame(
        {
            "n_muts": n_muts,
            "avg_vaf": vaf,
            "avg_n_mutant": n_mutant,
            "avg_n_total": n_total,
        },
        index=adata.obs_names,
    )


def get_vaf(adata):
    V = adata.to_df(layer="mutant")
    T = adata.to_df(layer="total")
    return V / T


def get_vaf_with_list(adata, alist):
    x = np.intersect1d(alist, adata.var_names)
    y = np.setdiff1d(alist, adata.var_names)
    df = get_vaf(adata[:, x])
    df[y] = None
    return df.transpose()


def bulk_total_mutant_vaf(adata):
    V = adata.to_df(layer="mutant")
    T = adata.to_df(layer="total")
    s1 = T.sum(axis=0).rename("total")
    s2 = V.sum(axis=0).rename("mutant")
    s3 = (s2 / s1).rename("vaf")
    return pd.concat([s1, s2, s3], axis=1)


def filter_mut_bulk_mutant_and_vaf_greater_than(adata, min_vaf=0.2, min_mutant=20):
    df = bulk_total_mutant_vaf(adata)
    good_muts = df[(df["vaf"] >= min_vaf) & (df["mutant"] >= min_mutant)].index
    keep_mut_by_list(adata, good_muts)


def filter_mut_total_coverage_must_be_greater_than(adata, min_total=7, min_cells=1):
    T = adata.layers["total"]
    good_muts = (T >= min_total).sum(axis=0) >= min_cells
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_vaf_greater_than_coverage_mutant_greater_than(
    adata, min_vaf=0.2, min_coverage_mutant=20, min_cells=1
):
    V = adata.to_df(layer="mutant")
    T = adata.to_df(layer="total")

    # good_muts_1 = (V >= min_coverage_mutant).sum(axis=0) >= min_cells
    df = V / T
    good_muts = ((df >= min_vaf) & (V >= min_coverage_mutant)).sum(axis=0) >= min_cells
    # a = adata.var_names.to_numpy()[good_muts_1]
    # b = adata.var_names.to_numpy()[good_muts_2]
    # good_muts = np.intersect1d(a, b)
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_coverage_mutant_must_be_greater_than(
    adata, min_coverage_mutant=3, min_cells=1
):
    V = adata.layers["mutant"]
    good_muts = (V >= min_coverage_mutant).sum(axis=0) >= min_cells
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def _pseudo_bulk_caller(sub_adata):
    genotype = 2 * np.ones(sub_adata.shape[1])
    mutant = sub_adata.layers["mutant"].sum(axis=0)
    total = sub_adata.layers["total"].sum(axis=0)

    genotype[total == mutant] = 3
    genotype[total > mutant] = 1
    genotype[mutant == 0] = 0
    genotype[total == 0] = 2

    return genotype, mutant, total


def _pseudo_bulk_caller_better(sub_adata, min_vaf, min_cell):
    genotype_final = 2 * np.ones(sub_adata.shape[1])
    genotype = (
        (sub_adata.layers["genotype"] == 1) | (sub_adata.layers["genotype"] == 3)
    ).sum(axis=0)
    mutant = sub_adata.layers["mutant"].sum(axis=0)
    total = sub_adata.layers["total"].sum(axis=0)
    vaf = mutant / total

    genotype_final[vaf > 0] = 0
    genotype_final[(vaf >= min_vaf) & (genotype >= min_cell)] = 1

    return genotype_final, mutant, total


def merge_cells_using(adata, using, min_vaf=0.4, min_cell=2):
    genotypes = []
    mutants = []
    totals = []
    indices = []
    for index, subgroup in adata.obs.groupby(using):
        genotype, mutant, total = _pseudo_bulk_caller_better(
            adata[subgroup.index, :], min_vaf, min_cell
        )
        indices.append(index)
        genotypes.append(genotype)
        mutants.append(mutant)
        totals.append(total)

    adata_merged = ad.AnnData(
        X=3 * np.ones((len(indices), adata.shape[1])),
        obs=pd.DataFrame.from_dict({"cells": indices}).set_index("cells"),
        var=adata.var,
        layers={
            "total": np.array(totals),
            "mutant": np.array(mutants),
            "genotype": np.array(genotypes),
        },
    )
    return adata_merged


def local_cluster_cells_then_merge_muts_pseudo_bulk(
    adata, by="mut", n_clusters=100, min_n_cells=5, attr="group"
):
    if by == "mut":
        clusters = tsc.ul.hclustering(adata.to_df())
    elif by == "cna":
        clusters = tsc.ul.hclustering(adata.obsm["cna"], metric="cosine")
    else:
        raise ValueError("Wrong `by` choice!")

    cluster = clusters[n_clusters]
    count = cluster.value_counts()
    cluster = cluster[cluster.isin(count[count >= min_n_cells].index)]
    cluster = cluster.apply(lambda x: f"G{x}")

    adata = adata[cluster.index, :]
    adata.obs[attr] = cluster

    genotypes = []
    mutants = []
    totals = []
    indices = []
    for index, subgroup in adata.obs.groupby(adata.obs[attr]):
        # genotype, mutant, total = _pseudo_bulk_caller_better(
        #     adata[subgroup.index, :], min_vaf=0.2, min_cell=1
        # )
        genotype, mutant, total = _pseudo_bulk_caller(adata[subgroup.index, :])

        indices.append(index)
        genotypes.append(genotype)
        mutants.append(mutant)
        totals.append(total)

    adata_merged = ad.AnnData(
        X=3 * np.ones((len(indices), adata.shape[1])),
        obs=pd.DataFrame.from_dict({"cells": indices}).set_index("cells"),
        var=adata.var,
        layers={
            "total": np.array(totals),
            "mutant": np.array(mutants),
            "genotype": np.array(genotypes),
        },
    )
    return adata_merged, adata


def filter_mut_mutant_must_present_in_at_least(adata, min_cells=2):
    """Remove mutations based on the mutant status.

    This function removes mutations that don't have the status of mutant
    in at least `min_cells` cells.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        The input readcount data.
    min_cells : :obj:`int`, optional
        Minimum number of cells, by default 1
    """

    G = adata.layers["genotype"]
    good_muts = ((G == 1) | (G == 3)).sum(axis=0) >= min_cells
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_mutant_in_groups_must_present_in_at_least(
    adata, min_groups=2, group_name="group", min_cells_within_group=1
):
    all_good_muts = []
    for index, subgroup in adata.obs.groupby(group_name):
        G = adata[subgroup.index, :].layers["genotype"]
        good_muts = ((G == 1) | (G == 3)).sum(axis=0) >= min_cells_within_group
        all_good_muts.append(good_muts)
    all_good_muts = np.array(all_good_muts)
    good_muts = all_good_muts.sum(axis=0) >= min_groups
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_reference_must_present_in_at_least(adata, min_cells=1):
    """Remove mutations based on the wild-type status.

    This function removes mutations that don't have the status of wild-type
    in at least `min_cells` cells. Note that mutations that are present in all cells
    will not be filtered out by this function.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        The input readcount data.
    min_cells : :obj:`int`, optional
        Minimum number of cells, by default 1
    """

    G = adata.layers["genotype"]
    good_muts = ((G == 0).sum(axis=0) >= min_cells) | (
        ((G == 1) | (G == 3)).sum(axis=0) == adata.shape[0]
    )
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_reference_in_groups_must_present_in_at_least(
    adata, min_groups=2, group_name="group", min_cells_within_group=1
):
    all_good_muts = []
    for index, subgroup in adata.obs.groupby(group_name):
        G = adata[subgroup.index, :].layers["genotype"]
        good_muts = ((G == 0).sum(axis=0) >= min_cells_within_group) & (
            ((G == 1) | (G == 3)).sum(axis=0) == 0
        )
        all_good_muts.append(good_muts)
    all_good_muts = np.array(all_good_muts)
    good_muts = all_good_muts.sum(axis=0) >= min_groups
    keep_mut_by_list(adata, adata.var_names.to_numpy()[good_muts])


def filter_mut_by_confidence_score(adata, min_p_value):
    T = adata.layers["total"]
    G = adata.layers["genotype"]

    # x = 0
    mut_to_remove = []
    for j, mut in enumerate(adata.var_names):
        mask = (G[:, j] == 1) | (G[:, j] == 3)
        g1 = T[:, j][mask]
        g2 = T[:, j][~mask]
        p_value = stats.ranksums(g1, g2)[1]
        if p_value < min_p_value:
            # x += 1
            mut_to_remove.append(mut)
            # tsc.logg.info(p_value, g1, g2)
        # if x == 10:
        #     break
    remove_mut_by_list(adata, mut_to_remove)


def remove_mut_by_list(adata, alist):
    """Remove a list of mutations from the data.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        The input readcount data.
    alist : :obj:`list`
        A list of mutations.
    """

    adata._inplace_subset_var(np.setdiff1d(adata.var_names, alist))
    tsc.logg.info(f"Matrix with n_obs × n_vars = {adata.shape[0]} × {adata.shape[1]}")


def keep_mut_by_list(adata, alist):
    adata._inplace_subset_var(np.intersect1d(alist, adata.var_names))
    tsc.logg.info(f"Matrix with n_obs × n_vars = {adata.shape[0]} × {adata.shape[1]}")


def remove_cell_by_list(adata, alist):
    """Remove a list of cells from the data.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        The input readcount data.
    alist : :obj:`list`
        A list of cells.
    """

    adata._inplace_subset_obs(np.setdiff1d(adata.obs_names, alist))
    tsc.logg.info(f"Matrix with n_obs × n_vars = {adata.shape[0]} × {adata.shape[1]}")


def keep_cell_by_list(adata, alist):
    adata._inplace_subset_obs(np.intersect1d(alist, adata.obs_names))
    tsc.logg.info(f"Matrix with n_obs × n_vars = {adata.shape[0]} × {adata.shape[1]}")


def get_germline_variants(adata, normal_cells, min_vaf=0.1, min_coverage=10):
    if not isinstance(normal_cells, list):
        raise "normal_cells must be a list"

    V = adata[normal_cells, :].to_df(layer="mutant")
    T = adata[normal_cells, :].to_df(layer="total")
    G = adata[normal_cells, :].to_df(layer="genotype")
    df = V / T
    good_muts = ((df >= min_vaf) & (V >= min_coverage) & ((G == 1) | (G == 3))).sum(
        axis=0
    ) >= 1
    return adata.var_names.to_numpy()[good_muts]


def build_scmatrix(adata):
    G = adata.layers["genotype"]
    adata.X[G == 0] = 0
    adata.X[(G == 1) | (G == 3)] = 1
    adata.X[G == 2] = 3
    adata.X = adata.X.astype("int")

    # M = adata.layers["mutant"]
    # T = adata.layers["total"]
    # adata.X[T != -1] = 3
    # adata.X[T > 0] = 0
    # adata.X[M > 0] = 1
    # adata.X = adata.X.astype("int")


def to_dendro(adata, filepath):
    Z = adata.layers["genotype"]
    Z[(Z == 1) | (Z == 3)] = 1
    Z[Z == 2] = 0
    X = adata.layers["mutant"]
    N = adata.layers["total"]
    pd.DataFrame(Z).to_csv(filepath + "/Z.csv")
    pd.DataFrame(X).to_csv(filepath + "/X.csv")
    pd.DataFrame(N).to_csv(filepath + "/N.csv")


def calc_relative_expression(adata, threshold):
    tpm = adata.to_df(layer="tpm")
    e_ij = np.log2(tpm / 10 + 1)

    # bad_cells = (tpm > 0).sum(axis=1) < 3000
    # tpm = tpm.loc[~bad_cells,:]
    # print(f'{sum(bad_cells)} cells filtered based on n_genes less than 3000')

    # bad_cells = e_ij[genes['housekeeping']].mean(axis=1) < 2.5
    # tpm = tpm.loc[~bad_cells,:]
    # e_ij = e_ij.loc[~bad_cells,:]
    # print(f'{sum(bad_cells)} cells filtered based on housekeeping less than 2.5')

    # bad_cells = e_ij[genes['mitochondria']].mean(axis=1) > 2.5
    # tpm = tpm.loc[~bad_cells,:]
    # e_ij = e_ij.loc[~bad_cells,:]
    # print(f'{sum(bad_cells)} cells filtered based on mitochondria greater than 2.5')

    bad_genes = (np.log2(tpm.mean(axis=0) + 1) < threshold).values
    e_ij = e_ij.loc[:, ~bad_genes]
    adata._inplace_subset_var(adata.var_names.to_numpy()[~bad_genes])
    tsc.logg.info(
        f"{sum(bad_genes)} genes filtered based on average accross cells less than"
        f" {threshold}"
    )

    tsc.logg.info(f"{adata.shape[0]} cells and {adata.shape[1]} genes retained")
    re_ij = e_ij - e_ij.groupby(lambda x: x.split("_")[0]).transform("mean")
    adata.layers["relative"] = re_ij.values


def add_vaf_to_layers(adata):
    V = adata.to_df(layer="mutant")
    T = adata.to_df(layer="total")
    df = V / T
    adata.layers["vaf"] = df


def statistics(adata):
    G = adata.layers["genotype"]
    t = adata.shape[0] * adata.shape[1]
    a = (G == 0).sum().sum()
    b = (G == 1).sum().sum()
    c = (G == 2).sum().sum()
    d = (G == 3).sum().sum()
    tsc.logg.info(f"size = {adata.shape[0]} × {adata.shape[1]}")
    tsc.logg.info(f"    HOM_REF = {a:6d} ({100*a/t:2.1f}%)")
    tsc.logg.info(f"    HET     = {b:6d} ({100*b/t:2.1f}%)")
    tsc.logg.info(f"    HOM_ALT = {d:6d} ({100*d/t:2.1f}%)")
    tsc.logg.info(f"    UNKNOWN = {c:6d} ({100*c/t:2.1f}%)")


def filter_private_mutations(df):
    df.drop(df.columns[df.sum() == 1], inplace=True)
