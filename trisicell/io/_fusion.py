import glob

import anndata as ad
import pandas as pd
from natsort import natsorted


def read_defuse(dirpath, min_span_count=10):
    data = []
    for file in natsorted(glob.glob(f"{dirpath}/*.tsv")):
        df = pd.read_table(file)
        name = file.split("/")[-1].replace(".tsv", "")
        df["fusion_id"] = df["gene_name1"] + "_" + df["gene_name2"]
        # df["fusion_id"] = (
        #     df["gene1"]
        #     + "_"
        #     + df["gene2"]
        #     + "_"
        #     + df["genomic_break_pos1"].astype(str)
        #     + "_"
        #     + df["genomic_break_pos2"].astype(str)
        # )
        df["cell_name"] = name
        df = df[df["span_count"] >= min_span_count]
        data.append(df)
    df = pd.concat(data)

    df = (
        df.pivot_table(index="cell_name", columns="fusion_id", values="span_count")
        .fillna(0)
        .astype(int)
    )
    df = df.loc[natsorted(df.index)]

    adata = ad.AnnData(1 * (df > 0), dtype=int)
    adata.layers["span_count"] = df.values
    return adata
