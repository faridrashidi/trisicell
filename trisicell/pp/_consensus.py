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
