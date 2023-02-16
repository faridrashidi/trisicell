import mudata as md
import scanpy as sc

import trisicell as tsc


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

    url = "https://ndownloader.figshare.com/files/39276362"
    sc.readwrite._check_datafile_present_and_download(
        "data/treated_igg_sw.h5md.gz", backup_url=url
    )
    mdata = md.read_h5mu("data/treated_igg_sw.h5md.gz")
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

    url = "https://ndownloader.figshare.com/files/22105673"
    sc.readwrite._check_datafile_present_and_download(
        "data/treated_igg_ss2.h5md.gz", backup_url=url
    )
    mdata = md.read_h5mu("data/treated_igg_ss2.h5md.gz")
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

    url = "https://ndownloader.figshare.com/files/22105670"
    sc.readwrite._check_datafile_present_and_download(
        "data/treated_actla4.h5md.gz", backup_url=url
    )
    mdata = md.read_h5mu("data/treated_actla4.h5md.gz")
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

    url = "https://ndownloader.figshare.com/files/22105682"
    sc.readwrite._check_datafile_present_and_download(
        "data/sublines_scrnaseq.h5md.gz", backup_url=url
    )
    mdata = md.read_h5mu("data/sublines_scrnaseq.h5md.gz")
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

    url = "https://ndownloader.figshare.com/files/22105691"
    sc.readwrite._check_datafile_present_and_download(
        "data/sublines_bwes.h5ad.gz", backup_url=url
    )
    adata = sc.read_h5ad("data/sublines_bwes.h5ad.gz")
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

    url = "https://ndownloader.figshare.com/files/22105685"
    sc.readwrite._check_datafile_present_and_download(
        "data/sublines_bwts.h5md.gz", backup_url=url
    )
    mdata = md.read_h5mu("data/sublines_bwts.h5md.gz")
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
