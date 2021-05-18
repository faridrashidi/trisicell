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
