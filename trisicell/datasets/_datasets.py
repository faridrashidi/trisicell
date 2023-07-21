import mudata as md
import scanpy as sc

import trisicell as tsc

url = "https://github.com/faridrashidi/trisicell/releases/download/v0.0.0"


def treated_igg_sw():
    """Trisicell treated mice (igg, seq-well) scRNAseq data.

    The size is n_cells × n_muts = 163 × 1453

    Returns
    -------
    :class:`mudata.MuData`

        A mudata with two modalities (`.mod`)

    Examples
    --------
    >>> mdata = tsc.datasets.treated_igg_sw()
    >>> mdata
    MuData object with n_obs × n_vars = 163 × 56854
    2 modalities
      mutation:	163 x 1453
        obs:	'group', 'Mouse_ID', 'Batch', 'Axl', 'Erbb3', 'Mitf', 'H2_K1', ...
        var:	'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', ...
        uns:	'trisicell_input', 'trisicell_output'
        layers:	'genotype', 'mutant', 'total'
      expression:	163 x 55401
        obs:	'uniquely_mapped_percent', 'num_splices', 'num_GCAG_splices', ...
        layers:	'fpkm', 'tpm'

    See Also
    --------
    :func:`trisicell.datasets.treated_actla4`.
    :func:`trisicell.datasets.treated_igg_ss2`.
    """

    name = "treated_igg_sw.h5md.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    mdata = md.read_h5mu(f"data/{name}")
    return mdata


def treated_igg_ss2():
    """Trisicell treated mice (igg, smart-seq2) scRNAseq data.

    The size is n_cells × n_muts = 163 × 1453

    Returns
    -------
    :class:`mudata.MuData`

        A mudata with two modalities (`.mod`)

    Examples
    --------
    >>> mdata = tsc.datasets.treated_igg_ss2()
    >>> mdata
    MuData object with n_obs × n_vars = 163 × 56854
    2 modalities
      mutation:	163 x 1453
        obs:	'group', 'Mouse_ID', 'Batch', 'Axl', 'Erbb3', 'Mitf', 'H2_K1', ...
        var:	'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', ...
        uns:	'trisicell_input', 'trisicell_output'
        layers:	'genotype', 'mutant', 'total'
      expression:	163 x 55401
        obs:	'uniquely_mapped_percent', 'num_splices', 'num_GCAG_splices', ...
        layers:	'fpkm', 'tpm'

    See Also
    --------
    :func:`trisicell.datasets.treated_actla4`.
    :func:`trisicell.datasets.treated_igg_sw`.
    """

    name = "treated_igg_ss2.h5md.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    mdata = md.read_h5mu(f"data/{name}")
    return mdata


def treated_actla4():
    """Trisicell treated mice (anti-ctla-4) scRNAseq data.

    The size is n_cells × n_muts = 508 × 3309

    Returns
    -------
    :class:`mudata.MuData`

        A mudata with two modalities (`.mod`)

    Examples
    --------
    >>> mdata = tsc.datasets.treated_actla4()
    >>> mdata
    MuData object with n_obs × n_vars = 508 × 58710
    2 modalities
      mutation:	508 x 3309
        obs:	'group', 'Mouse_ID', 'Batch', 'Axl', 'Erbb3', 'Mitf', 'H2-K1', ...
        var:	'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', ...
        uns:	'trisicell_input', 'trisicell_output'
        layers:	'genotype', 'mutant', 'total'
      expression:	508 x 55401
        obs:	'uniquely_mapped_percent', 'num_splices', 'num_GCAG_splices', ...
        layers:	'fpkm', 'tpm'

    See Also
    --------
    :func:`trisicell.datasets.treated_igg_ss2`.
    :func:`trisicell.datasets.treated_igg_sw`.
    """

    name = "treated_actla4.h5md.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    mdata = md.read_h5mu(f"data/{name}")
    return mdata


def sublines_scrnaseq():
    """Trisicell sublines scRNAseq data.

    The size is n_cells × n_muts = 175 × 450

    Returns
    -------
    :class:`mudata.MuData`

        A mudata with two modalities (`.mod`)

    Examples
    --------
    >>> mdata = tsc.datasets.sublines_scrnaseq()
    >>> mdata
    MuData object with n_obs × n_vars = 175 × 55851
    2 modalities
      expression: 175 x 55401
        obs:	'cells', 'uniquely_mapped_percent', 'num_splices', ...
        layers:	'fpkm', 'tpm'
      mutation: 175 x 450
        obs:	'cells', 'clone', 'group', 'group_color', 'is_red', 'is_sub', ...
        var:	'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', ...
        layers:	'genotype', 'mutant', 'total', 'trisicell_input', 'trisicell_output'

    See Also
    --------
    :func:`trisicell.datasets.sublines_bwes`.
    :func:`trisicell.datasets.sublines_bwts`.
    """

    name = "sublines_scrnaseq.h5md.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    mdata = md.read_h5mu(f"data/{name}")
    return mdata


def sublines_bwes():
    """Trisicell sublines bWES data.

    The size is n_sublines × n_muts = 24 × 6653

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.var` contains information about the mutations.

            - `.layers['trisicell_input']` the binary input genotype matrix used as
                input to the Trisicell.
            - `.layers['trisicell_output']` the binary input genotype matrix inferred by
                Trisicell-boost(SCITE).
            - `.layers['genotype']` noisy genotype matrix, 0: reference, 1: heterozygous
                2: unknown and 3: homozygous_alt.
            - `.layers['mutant']` number of mutant reads.
            - `.layers['total']` number of total reads.

    Examples
    --------
    >>> adata = tsc.datasets.sublines_bwes()
    >>> adata
    AnnData object with n_obs × n_vars = 24 × 6653
        var: 'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', 'position', ...
        layers: 'genotype', 'mutant', 'total', 'trisicell_input', 'trisicell_output'

    See Also
    --------
    :func:`trisicell.datasets.sublines_scrnaseq`.
    :func:`trisicell.datasets.sublines_bwts`.
    """

    name = "sublines_bwes.h5ad.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    adata = sc.read_h5ad(f"data/{name}")
    return adata


def sublines_bwts():
    """Trisicell sublines bWTS data.

    The size is n_cells × n_muts = 33 × 536

    Returns
    -------
    :class:`mudata.MuData`

        A mudata with two modalities (`.mod`)

    Examples
    --------
    >>> mdata = tsc.datasets.sublines_bwts()
    >>> mdata
    MuData object with n_obs × n_vars = 33 × 55851
    2 modalities
      expression: 175 x 55401
        obs:	'cells', 'uniquely_mapped_percent', 'num_splices', ...
        layers:	'fpkm', 'tpm'
      mutation:	33 x 536
        obs:	'cells', 'clone', 'mps', 'zscore', 'group', 'Axl', 'Mitf', 'day', ...
        var:	'kind', 'amino_acid_change', 'ensemble', 'gene', 'chrom', 'position',...
        layers:	'genotype', 'mutant', 'total', 'trisicell_input', 'trisicell_output'

    See Also
    --------
    :func:`trisicell.datasets.sublines_bwes`.
    :func:`trisicell.datasets.sublines_scrnaseq`.
    """

    name = "sublines_bwts.h5md.gz"
    sc.readwrite._check_datafile_present_and_download(
        f"data/{name}", backup_url=f"{url}/{name}"
    )
    mdata = md.read_h5mu(f"data/{name}")
    return mdata


def example(is_expression=False):
    """Return an example for sanity checking and playing with Trisicell.

    is_expression : :obj:`bool`, optional
        Returns the expression dataset instead of the genotype one, by default False

    Returns
    -------
    :class:`anndata.AnnData`
        An object that cells are in `.obs` and mutations are in `.var`.
    """

    if is_expression:
        return tsc.io.read(
            tsc.ul.get_file("trisicell.datasets/data/expression.h5ad.gz")
        )
    else:
        return tsc.io.read(tsc.ul.get_file("trisicell.datasets/data/genotype.h5ad.gz"))


def test():
    df = tsc.io.read(tsc.ul.get_file("trisicell.datasets/test/test.tsv"))
    return df


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
            - `.layers['solution_fig7a']` is the solution presented in Figure 7a of
                PhISCS paper.
            - `.layers['solution_fig7b']` is the solution presented in Figure 7b of
                PhISCS paper.
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
    return adata


def high_grade_serous_ovarian_cancer_3celllines():
    """High Grade Serous Ovarian Cancer (3 cell lines).

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

            - `.obs['clone_id']` is the clone id to which the cell is assigned in
                Figure 3H.
            - `.obs['group_color']` is unique colors for each 'clone_id'.
            - `.obs['cell_name']` is a new name for each cell based on the
                'group_color'.
            - `.layers['mutant']` is the number of mutant reads at each locus in each
                cell.
            - `.layers['total']` is the total number of reads at each locus in each
                cell.
            - `.layers['ground']` is the solution inferred in Figure 3H of the original
                paper.
            - `.uns['params_ground']` is parameters inferred by comparing ground and
                noisy matrices.
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
