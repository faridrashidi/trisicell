import trisicell as tsc


def binarym_filter_private_mutations(df):
    df.drop(df.columns[df.sum() == 1], axis=1, inplace=True)


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
