import pandas as pd

import trisicell as tsc


def example(kind):
    """A toy example for sanity checking and playing with Trisicell.

    Parameters
    ----------
    kind : :obj:`str`
        The options are:

            - `'genotype'`: This option will return information about genotype matrix.
            - `'expression'`: This option will return information about expression matrix.

    Returns
    -------
    :class:`anndata.AnnData`
        An object that cells are in `.obs` and mutations are in `.var`.

    Raises
    ------
    ValueError
        If `kind` is not `genotype` or `expression`.
    """

    if kind.lower() == "genotype":
        return tsc.io.read(tsc.ul.get_file("trisicell.datasets/data/genotype.h5ad.gz"))
    elif kind.lower() == "expression":
        return tsc.io.read(
            tsc.ul.get_file("trisicell.datasets/data/expression.h5ad.gz")
        )
    else:
        raise ValueError("kind must be either genotype or expression")


def test():
    df = tsc.io.read(tsc.ul.get_file("trisicell.datasets/test/test.tsv"))
    return df


def melanoma20():
    """Mouse Melanoma dataset with 20 sublines.

    This dataset was introduced in :cite:`Wolf_2019` and was used in:

    * :cite:`PhISCS-BnB` Figure 1.

    The size is n_cells × n_muts = 20 × 2367

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

        `.layers['solution_fig1']` is the solution presented in Figure 1 of PhISCS-BnB paper.
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/melanoma20.h5ad"))
    return adata


def colorectal1():
    """Human Colorectal Cancer (Patient 1).

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`B-SCITE` Figure 8a.
    * :cite:`SiFit` Figure 6.
    * :cite:`SPhyR` Table 1.
    * :cite:`SiCloneFit` Figure 3.

    The size is n_cells × n_muts = 72 × 12

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/colorectal1.rc.h5ad"))
    # FIXME: extract (SiFit 178 × 16).
    # https://github.com/algo-cancer/PhyloM/blob/master/Data/README.md
    pass


def colorectal2(readcount=False):
    """Human Colorectal Cancer (Patient 2).

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`PhISCS` Figure 7.
    * :cite:`B-SCITE` Figure 8b.
    * :cite:`SiCloneFit` Figure 4.
    * :cite:`SCARLET` Figure 4.

    The size is n_cells × n_muts = 78 × 25

    Parameters
    ----------
    readcount : :obj:`str`
        Return the readcount information of the data.

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.layers['solution_fig7a']` is the solution presented in Figure 7a of PhISCS paper.
            - `.layers['solution_fig7b']` is the solution presented in Figure 7b of PhISCS paper.
            - `.uns['params_fig7a']` is parameters used as input to get 'solution_7a'.
            - `.uns['params_fig7b']` is parameters used as input to get 'solution_7b'.
            - `.var` includes information of the bulk samples.
    """

    if readcount:
        adata = tsc.io.read(
            tsc.ul.get_file("trisicell.datasets/real/colorectal2.rc.h5ad")
        )
    else:
        adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/colorectal2.h5ad"))
    # FIXME: (86 x 25 in B-SCITE) (182 x 36 SiCloneFit)
    # https://github.com/cbg-ethz/infSCITE/blob/master/pat_2.csv
    return adata


def colorectal3():
    """Human Colorectal Cancer.

    This dataset was introduced in :cite:`Wu_2016` and was used in:

    * :cite:`SiFit` Figure 5.

    The size is n_cells × n_muts = 48 × 77

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/colorectal3.h5ad"))
    # TODO: extract.
    return adata


def acute_lymphocytic_leukemia1():
    """Human Acute Lymphocytic Leukemia dataset (Patient 1).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`B-SCITE` Figure 5.
    * :cite:`infSCITE` Figure S16.

    The size is n_cells × n_muts = 111 × 20

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia1.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia2():
    """Human Acute Lymphocytic Leukemia dataset (Patient 2).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`PhISCS` Figure 9.
    * :cite:`B-SCITE` in Figure 6.
    * :cite:`infSCITE` Figure S17.
    * :cite:`Phyolin` Table 2.

    The size is n_cells × n_muts = 102 × 16

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.layers['solution_fig9']` is the solution presented in Figure 9 of PhISCS paper.
            - `.uns['params_fig9']` is parameters used as input to get 'solution_fig9'.
            - `.var` includes information of the bulk samples.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia2.h5ad")
    )
    # FIXME: 115 x 16 in B-SCITE?
    return adata


def acute_lymphocytic_leukemia3():
    """Human Acute Lymphocytic Leukemia dataset (Patient 3).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S18.
    * :cite:`SCIPhI` Figure 5.
    * :cite:`ScisTree` Figure S3.

    The size is n_cells × n_muts = 150 × 49

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia3.h5ad")
    )
    # FIXME: 255 x 49 in SCIPhI and ScisTree?
    return adata


def acute_lymphocytic_leukemia4():
    """Human Acute Lymphocytic Leukemia dataset (Patient 4).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S19.
    * :cite:`gpps` Figure 3.
    * :cite:`SASC` Figure 6.

    The size is n_cells × n_muts = 143 × 78

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia4.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia5():
    """Human Acute Lymphocytic Leukemia dataset (Patient 5).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S20.
    * :cite:`SASC` Figure 7.

    The size is n_cells × n_muts = 96 × 105

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia5.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia6():
    """Human Acute Lymphocytic Leukemia dataset (Patient 6).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S21.
    * :cite:`Phyolin` Table 2.

    The size is n_cells × n_muts = 146 × 10

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/acute_lymphocytic_leukemia6.h5ad")
    )
    return adata


def tnbc():
    """Triple-negative Breast Cancer.

    This dataset was introduced in :cite:`Wang_2014` and was used in:

    * :cite:`SCIPhI` Figure 4.
    * :cite:`B-SCITE` Figure 7.
    * :cite:`TRaIT` Figure 6.

    The size is n_cells × n_muts = 16 × 20

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

    Examples
    --------
    >>> adata = tsc.datasets.tnbc()
    >>> df_in = adata.to_df()
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/tnbc.h5ad"))
    return adata


def erbc():
    """Oestrogen-receptor-positive (ER+) Breast Cancer.

    This dataset was introduced in :cite:`Wang_2014` and was used in:

    * :cite:`SCITE` Figure S8 and S9.
    * :cite:`infSCITE` Figure S15.
    * :cite:`gpps` Figure 1.

    The size is n_cells × n_muts = 47 × 40

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Examples
    --------
    >>> adata = tsc.datasets.erbc()
    >>> df_in = adata.to_df()
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/erbc.h5ad"))
    return adata


def muscle_invasive_bladder():
    """Muscle Invasive Bladder Cancer.

    This dataset was introduced in :cite:`Li_2012` and was used in:

    * :cite:`OncoNEM` Figure 6B.

    The size is n_cells × n_muts = 44 × 443

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_onconem']` is parameters inferred by OncoNEM.

    Examples
    --------
    >>> adata = tsc.datasets.muscle_invasive_bladder()
    >>> df_in = adata.to_df()
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/muscle_invasive_bladder.h5ad")
    )
    return adata


def renal_cell_carcinoma():
    """Clear-cell Renal-cell Carcinoma.

    This dataset was introduced in :cite:`Xu_2012` and was used in:

    * :cite:`SCITE` Figure S6 and S7.
    * :cite:`infSCITE` Figure S14.

    The size is n_cells × n_muts = 17 × 35

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Examples
    --------
    >>> adata = tsc.datasets.renal_cell_carcinoma()
    >>> df_in = adata.to_df()
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/renal_cell_carcinoma.h5ad")
    )
    return adata


def myeloproliferative_neoplasms712():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`OncoNEM` Figure 6D.

    The size is n_cells × n_muts = 58 × 712

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_onconem']` is parameters inferred by OncoNEM.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/myeloproliferative_neoplasms712.h5ad")
    )
    return adata


def myeloproliferative_neoplasms78():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`SCITE` Figure S5.

    The size is n_cells × n_muts = 58 × 78

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

    Notes
    -----
    The original dataset contains 712 mutations but 78 ones were considered as
    non-synonymous mutations from the full data.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/myeloproliferative_neoplasms78.h5ad")
    )
    return adata


def myeloproliferative_neoplasms18():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`SCITE` Figure S2, S3 and S4.
    * :cite:`Kim_2014` Figure 1.
    * :cite:`infSCITE` Figure S13.
    * :cite:`gpps` Figure 2.

    The size is n_cells × n_muts = 58 × 18

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Notes
    -----
    The original dataset contains 712 mutations but 18 ones were considered as
    cancer related mutations from the full data.
    """

    adata = tsc.io.read(
        tsc.ul.get_file("trisicell.datasets/real/myeloproliferative_neoplasms18.h5ad")
    )
    return adata


def high_grade_serous_ovarian_cancer1():
    """Triple-negative Breast Cancer (Patient 2).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S22.

    The size is n_cells × n_muts = 588 × 37

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file(
            "trisicell.datasets/real/high_grade_serous_ovarian_cancer1.h5ad"
        )
    )
    # TODO: extract.
    return adata


def high_grade_serous_ovarian_cancer2():
    """Triple-negative Breast Cancer (Patient 3).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S23.

    The size is n_cells × n_muts = 672 × 60

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file(
            "trisicell.datasets/real/high_grade_serous_ovarian_cancer2.h5ad"
        )
    )
    # TODO: extract.
    return adata


def high_grade_serous_ovarian_cancer3():
    """Triple-negative Breast Cancer (Patient 9).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S24.

    The size is n_cells × n_muts = 420 × 37

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    adata = tsc.io.read(
        tsc.ul.get_file(
            "trisicell.datasets/real/high_grade_serous_ovarian_cancer3.h5ad"
        )
    )
    # TODO: extract.
    return adata


def high_grade_serous_ovarian_cancer_3celllines():
    """Triple-negative Breast Cancer ().

    This dataset was introduced in :cite:`Laks_2019`
    and preprocessed in :cite:`McPherson_2019`.
    The phylogeny is presented in Figure 3H.

    The size is n_cells × n_muts = 891 × 13666

    Note that 402 mutations were deleted during evolution
    in the inferred tree by original study So, here they were filtered out
    (meaning they are 14068 mutations in total in the original study).

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy
        (it is obtained based on the number of mutant/total reads).

            - `.obs['clone_id']` is the clone id to which the cell is assigned in Figure 3H.
            - `.obs['group_color']` is unique colors for each 'clone_id'.
            - `.obs['cell_name']` is a new name for each cell based on the 'group_color'.
            - `.layers['mutant']` is the number of mutant reads at each locus in each cell.
            - `.layers['total']` is the total number of reads at each locus in each cell.
            - `.layers['ground']` is the solution inferred in Figure 3H of the original paper.
            - `.uns['params_ground']` is parameters inferred by comparing ground and noisy matrices.
            - `.var` includes information of the bulk samples.

    Examples
    --------
    >>> sc = tsc.datasets.high_grade_serous_ovarian_cancer_3celllines()
    >>> print(sc)
    AnnData object with n_obs × n_vars = 891 × 13666
        obs: 'clone_id', 'group_color', 'cell_name'
        uns: 'params_ground'
        layers: 'ground', 'mutant', 'total'
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/ovarian.h5ad.gz"))
    return adata
