import glob

import anndata as ad
import pandas as pd
from natsort import natsorted


def read_defuse(dirpath):
    data = []
    raw = []
    for file in natsorted(glob.glob(f"{dirpath}/*.tsv")):
        df = pd.read_table(file)
        if df.shape[0] == 0:
            continue
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
        # df["fusion_id"] = (
        #     df["gene_name1"]
        #     + "_"
        #     + df["gene_name2"]
        #     + "_"
        #     + df["genomic_break_pos1"].astype(str)
        #     + "_"
        #     + df["genomic_break_pos2"].astype(str)
        # )
        df["cell_name"] = name
        raw.append(df)
        df = df[df["adjacent"] == "N"]
        df["gene_chromosome1"] = df["gene_chromosome1"].astype(str)
        df["gene_chromosome2"] = df["gene_chromosome2"].astype(str)
        df = df[df["gene_chromosome1"] != df["gene_chromosome2"]]
        data.append(df)
    data = pd.concat(data)
    raw = pd.concat(raw)
    raw

    data = (
        data.pivot_table(index="cell_name", columns="fusion_id", values="span_count")
        .fillna(0)
        .astype(int)
    )
    data = data.loc[natsorted(data.index)]

    adata = ad.AnnData(1 * (data > 0), dtype=int)
    adata.layers["span_count"] = data.values
    return adata
