import click

import trisicell as tsc


@click.command(short_help="Calculate scores between two trees.")
@click.argument(
    "ground_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "inferred_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def score(ground_file, inferred_file):
    """Scores of Anscenstor-descendant, Different-lineage and MLTD.

    trisicell score ground.CFMatrix inferred.CFMatrix
    """

    tsc.settings.verbosity = "info"

    df_g = tsc.io.read(ground_file)
    df_s = tsc.io.read(inferred_file)

    gs = tsc.tl.gs(df_g, df_s)
    ad = tsc.tl.ad(df_g, df_s)
    dl = tsc.tl.dl(df_g, df_s)
    mltd = tsc.tl.mltd(df_g, df_s)["normalized_similarity"]
    tpted = tsc.tl.tpted(df_g, df_s)
    caset = tsc.tl.caset(df_g, df_s)
    disc = tsc.tl.disc(df_g, df_s)
    mp3 = tsc.tl.mp3(df_g, df_s)
    rf = tsc.tl.rf(df_g, df_s)

    tsc.logg.info(
        f"GS={gs:0.4f}\nAD={ad:0.4f}\nDL={dl:0.4f}\n"
        f"MLTSM={mltd:0.4f}\nTPTED={tpted:0.4f}\n"
        f"CASet={caset:0.4f}\nDISC={disc:0.4f}\nMP3={mp3:0.4f}\n"
        f"RF={rf:0.4f}"
    )

    return None
