import numpy as np
import pandas as pd
import seaborn as sns

import trisicell as tsc


def difference_between(dna, rna):
    dna_temp = (
        dna.var["chrom"].astype(str)
        + ":"
        + dna.var["position"].astype(str)
        + ":"
        + dna.var["reference"].astype(str)
        + ":"
        + dna.var["alteration"].astype(str)
    )
    rna_temp = (
        rna.var["chrom"].astype(str)
        + ":"
        + rna.var["position"].astype(str)
        + ":"
        + rna.var["reference"].astype(str)
        + ":"
        + rna.var["alteration"].astype(str)
    )
    a = np.intersect1d(dna_temp.values, rna_temp.values)
    b = np.setdiff1d(dna_temp.values, rna_temp.values)
    c = np.setdiff1d(rna_temp.values, dna_temp.values)
    tsc.logg.info(f"{len(b)}, inter={len(a)}, {len(c)}")

    b = dna.to_df(layer="mutant").T.loc[dna_temp[dna_temp.isin(b)].index]
    c = rna.to_df(layer="mutant").T.loc[rna_temp[rna_temp.isin(c)].index]

    return dna_temp[dna_temp.isin(a)].index, b, c


def flips_in_cells(adata, df_in, df_out):
    c01 = ((df_in == 0) & (df_out == 1)).sum(axis=1)
    c10 = ((df_in == 1) & (df_out == 0)).sum(axis=1)
    c31 = ((df_in == 3) & (df_out == 1)).sum(axis=1)
    c30 = ((df_in == 3) & (df_out == 0)).sum(axis=1)
    c01 = c01 / c01.sum()
    c10 = c10 / c10.sum()
    c31 = c31 / c31.sum()
    c30 = c30 / c30.sum()
    c01 = c01.rename("f_0_1")
    c10 = c10.rename("f_1_0")
    c31 = c31.rename("f_3_1")
    c30 = c30.rename("f_3_0")
    df = pd.DataFrame([c01, c10, c31, c30]).T
    df["f_0_1_color"] = "#E41A1C"
    df["f_1_0_color"] = "#FF7F00"
    df["f_3_1_color"] = "#FDBF6F"
    df["f_3_0_color"] = "#FB9A99"
    adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True)


def get_database(species):
    data = []
    temp = {}
    with open(f"/home/frashidi/database/{species}/ensembl/{species}.rsem") as fin:
        i = -2
        for line in fin:
            i += 1
            line = line.strip()
            if i == -1:
                continue
            if i % 6 == 0:
                temp = {}
                temp["trans_ens"] = line.split("\t")[0]
                temp["trans_sym"] = line.split("\t")[1]
            if i % 6 == 1:
                temp["gene_ens"] = line.split("\t")[0]
                temp["gene_sym"] = line.split("\t")[1]
            if i % 6 == 2:
                temp["chrom"] = line
            if i % 6 == 3:
                temp["strand"] = line.split(" ")[0]
            if i % 6 == 4:
                temp["n_exons"] = int(line.split(" ")[0])
                temp["regions"] = line.split(" ")[1:]
                temp["start"] = int(line.split(" ")[1])
                temp["end"] = int(line.split(" ")[-1])
            if i % 6 == 5:
                details = line.split("; ")
                for det in details:
                    temp[det.split(" ")[0]] = det.split(" ")[1]
                data.append(temp)
    database = pd.DataFrame(data)
    return database


def translate_gene(alist, kind):
    mm10 = get_database("mm10")
    mm10 = mm10.set_index("gene_sym")["gene_ens"].drop_duplicates().to_dict()
    hg19 = get_database("hg19")
    hg19 = hg19.set_index("gene_sym")["gene_ens"].drop_duplicates().to_dict()
    if kind == "mm10":
        translator = mm10
    else:
        translator = hg19
    result = []
    for d in alist:
        if d in translator:
            result.append(f"{translator[d]}_{d}")
        else:
            tsc.logg.info(d)
    return result


def split_mut(mut):
    try:
        ref = mut.split(".")[-2]
        pos = mut.split(".")[-3]
        chrom = mut.split(".")[-4]
        gene = mut.split(".chr")[0].split("_")[1]
        ens = mut.split(".chr")[0].split("_")[0]
        alt = mut.split(".")[-1]
        return ens, gene, chrom, pos, ref, alt
    except Exception:
        return None, None, None, None, None, None


def colorize_one(sery, color="#984EA3"):
    x = sery.sort_values()
    return pd.Series(
        sns.color_palette(f"light:{color}", len(x)).as_hex(), index=x.index
    )


def colorize_two(sery):
    x = sery.sort_values()
    a = (x <= 0).sum()
    b = (0 < x).sum()
    a = sns.color_palette("Reds_r", a).as_hex()
    b = sns.color_palette("Blues", b).as_hex()
    return pd.Series(a + b, index=x.index)
