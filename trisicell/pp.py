def filter(df):
    """Filter input.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Input dataframe.

    Returns
    -------
    :class:`pandas.DataFrame`
        Output dataframe.
    """
    return df[df > 0]
