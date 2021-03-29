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
    df = tsc.io.read(tsc.ul.get_file("trisicell.datasets/test/test.SC"))
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
    """Human Colorectal Cancer dataset patient 1.

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`B-SCITE` Figure 8a.

    The size is n_cells × n_muts = 72 × 12

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    # TODO: extract.
    pass


def colorectal2():
    """Human Colorectal Cancer dataset patient 2.

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`PhISCS` Figure 7.
    * :cite:`B-SCITE` Figure 8b.

    The size is n_cells × n_muts = 78 × 25

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

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/colorectal2.h5ad"))
    # TODO: 86 x 25 in B-SCITE?
    return adata


def acute_lymphocytic_leukemia1():
    """Human Acute Lymphocytic Leukemia dataset.

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`B-SCITE` Figure 5.

    The size is n_cells × n_muts = 111 × 20

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    # TODO: extract.
    pass


def acute_lymphocytic_leukemia2():
    """Human Acute Lymphocytic Leukemia dataset.

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`PhISCS` Figure 9.
    * :cite:`B-SCITE` in Figure 6.

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
    # TODO: 115 x 16 in B-SCITE?
    return adata


def acute_lymphocytic_leukemia3():
    # TODO: extract.
    pass


def acute_lymphocytic_leukemia4():
    # TODO: extract.
    pass


def acute_lymphocytic_leukemia5():
    # TODO: extract.
    pass


def acute_lymphocytic_leukemia6():
    # TODO: extract.
    pass


def myeloproliferative_neoplasms18():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`PhISCS` Figure 9.

    The size is n_cells × n_muts = 58 × 18

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.layers['solution_fig9']` is the solution presented in Figure 9 of PhISCS paper.
            - `.uns['params_fig9']` is parameters used as input to get 'solution_fig9'.
            - `.var` includes information of the bulk samples.
    """

    pass


def triple_negative_breast_cancer():
    """Triple-negative Breast Cancer.

    This dataset was introduced in []_ and was used in:
    * :cite:`B-SCITE` Figure 7.

    The size is n_cells × n_muts = 16 × 18

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """

    # TODO: extract.
    pass


def high_grade_serous_ovarian_cancer():
    """Triple-negative Breast Cancer.

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
    >>> sc = tsc.datasets.high_grade_serous_ovarian_cancer()
    >>> print(sc)
    AnnData object with n_obs × n_vars = 891 × 13666
        obs: 'clone_id', 'group_color', 'cell_name'
        uns: 'params_ground'
        layers: 'ground', 'mutant', 'total'
    """

    adata = tsc.io.read(tsc.ul.get_file("trisicell.datasets/real/ovarian.h5ad.gz"))
    return adata
