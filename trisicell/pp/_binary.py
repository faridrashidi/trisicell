import trisicell as tsc


def binarym_filter_private_mutations(df):
    df.drop(df.columns[df.sum() == 1], axis=1, inplace=True)


def binarym_filter_clonal_mutations(df):
    x = (df == 1).sum()
    x = x[x == df.shape[0]]
    df.drop(x.index, axis=1, inplace=True)


def binarym_filter_nonsense_mutations(df, alt_in=2, ref_in=1):
    df.drop(
        df.columns[
            ~(
                ((df == 1).sum() >= alt_in)
                & (
                    ((df == 0).sum() >= ref_in)
                    | ((df == 1).sum() == df.shape[0])
                    | ((df == 1).sum() >= (df == 3).sum())
                )
            )
        ],
        axis=1,
        inplace=True,
    )


def binarym_statistics(df):
    t = df.shape[0] * df.shape[1]
    a = (df == 0).sum().sum()
    b = (df == 1).sum().sum()
    d = (df == 3).sum().sum()
    tsc.logg.info(f"size = {df.shape[0]} × {df.shape[1]}")
    tsc.logg.info(f"    REF     = {a:6d} ({100*a/t:2.1f}%)")
    tsc.logg.info(f"    HET     = {b:6d} ({100*b/t:2.1f}%)")
    tsc.logg.info(f"    UNKNOWN = {d:6d} ({100*d/t:2.1f}%)")


def consensus_combine(df):
    """Combine cells in genotype matrix.

    This function combines the replicates or cells that have
    the same origin prior to running Trisicell-Cons. The replicates
    or cells that are supposed to be merged must be designated
    with `_`. For instance:

    input: {`{Cell1}_{ID1}`, `{Cell1}_{ID2}`, `{Cell2}_{ID1}`, `{Cell2}_{ID2}`}.

    output: {`{Cell1}`, `{Cell2}`}.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        The input genotype matrix in conflict-free format.

    Returns
    -------
    :class:`pandas.DataFrame`
        The combine genotype matrix.
    """

    df2 = df.groupby(df.index.str.split("_").str[0]).transform("prod")
    df2 = df2.groupby(df2.index.str.split("_").str[0]).first()
    return df2
