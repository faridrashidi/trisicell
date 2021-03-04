import pandas as pd


def read(path):
    """Read input data.

    Parameters
    ----------
    path : str
        The path to the input file.

    Returns
    -------
    :class:`pandas.DataFrame`
        The input dataframe.
    """
    df = pd.read_csv(path)
    return df
