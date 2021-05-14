import glob
import operator
import os
import pickle
import sys

import anndata as ad
import numpy as np
import pandas as pd
from cyvcf2 import VCF

import trisicell as tsc


def read(filepath):
    """Read genotype matrix and read-count matrix.

    The genotype matrix must be in the in format of :class:`pandas.DataFrame`
    The read-count matrix must be in the format of :class:`anndata.AnnData`.

    Parameters
    ----------
    filepath : :obj:`str`
        The path to the file. The extension must be one of
        [`.SC`, `.CFMatrix`, `.h5ad`, `.h5ad.gz`]

    Returns
    -------
    :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        Depends on the format of the input file the output type is different.

    Raises
    ------
    ValueError
        If the extension is not one of the above.
    """

    ext = os.path.splitext(filepath)[-1]
    if ext in [".SC", ".CFMatrix", ".before_FP_FN_NA", ".tsv"]:
        sc = pd.read_table(filepath, index_col=0)
        if len(sc.columns) != len(set(sc.columns)):
            raise RuntimeError("Mutation ids must be unique!")
        return sc
    elif ext in [".h5ad", ".gz"]:
        return ad.read(filepath)
    else:
        raise ValueError("Extension is wrong!")


def write(obj, filepath):
    """Write genotype matrix or read-count matrix into a file.

    Parameters
    ----------
    obj : :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        The input object which is going to be written in a file.
    filepath : :obj:`str`
        The file path where the `obj` must be written in.

    Raises
    ------
    ValueError
        If `obj` is not an instance of:
            :class:`pandas.DataFrame`
            :class:`anndata.AnnData`
    """

    if isinstance(obj, pd.DataFrame):
        obj.index.name = "cellIDxmutID"
        obj.to_csv(filepath, sep="\t")
    elif isinstance(obj, ad.AnnData):
        obj.write(filepath + ".h5ad.gz", compression="gzip")
    else:
        raise ValueError("Object instance is wrong!")


def read_gatk(
    folderpath,
    genome_id,
    only_snp=False,
    only_exonic=True,
    filter_private=False,
    filter_dbsnp=True,
    filter_indels_greater_than_2=True,
):
    """Read GATK genotype output to AnnData matrix.

    Parameters
    ----------
    folderpath : :obj:`str`
        The path to the folder contains following files:

        1) `jointcalls.g.vcf`
        2) `exonic.tsv`
        3) `readstat.tsv`
    genome_id : :obj:`str`
        One of `hg19`, `mm10`
    only_snp : :obj:`bool`, optional
        If only SNPs must be called (i.e. Indels are ignored), by default False
    only_exonic : :obj:`bool`, optional
        Only exonic mutations are considered, by default True
    filter_private : :obj:`bool`, optional
        Filter private mutation (i.e. only present in one cell), by default False
    filter_dbsnp : :obj:`bool`, optional
        Filter mutations that are present in dbSNP, by default True
    filter_indels_greater_than_2 : :obj:`bool`, optional
        Filter Indels that the length is greater than 2, by default True

    Returns
    -------
    :class:`anndata.AnnData`
        The AnnData object which includes layers of mutant, total and genotype
    """

    def get_germline(species):
        filename = f"/home/frashidi/database/{species}/dbsnp/{species}.pickle"
        obj = pickle.load(open(filename, "rb"))
        return obj

    def get_database(species):
        data = []
        with open(f"/home/frashidi/database/{species}/ensembl/{species}.rsem") as fin:
            i = -2
            for line in fin:
                i += 1
                line = line.strip()
                if i == -1:
                    continue
                if i % 6 == 1:
                    symbol = line.split("\t")[1]
                if i % 6 == 2:
                    chrom = line
                if i % 6 == 3:
                    strand = line.split(" ")[0]
                if i % 6 == 5:
                    details = line.split(";")
                    ensg = details[0].split("gene_id ")[1].replace('"', "")
                    enst = details[1].split("transcript_id ")[1].replace('"', "")
                    data.append(
                        {
                            "enst_symbol": f"{enst}_{symbol}",
                            "ensg_symbol": f"{ensg}_{symbol}",
                        }
                    )
        database = pd.DataFrame(data).set_index("enst_symbol")["ensg_symbol"].to_dict()
        return database

    def get_exonic(exonic_file, species):
        database = get_database(species)
        df = pd.read_table(exonic_file, header=None)
        df["muts"] = df.apply(lambda x: f"{x[1]}.{x[2]}.{x[3]}.{x[4]}", axis=1)
        df["gene"] = df[0].map(database)
        df["kind"] = df[5]
        df["amino_acid_change"] = df[6]
        df = df.drop([0, 1, 2, 3, 4, 5, 6], axis=1)
        df = df.drop_duplicates("muts")
        return df.set_index("muts")

    vcf_file = f"{folderpath}/jointcalls.g.vcf"
    exonic_file = f"{folderpath}/exonic.tsv"
    if genome_id == "human":
        species = "hg19"
        germline = get_germline(species)
    else:
        species = "mm10"
        germline = get_germline(species)

    if not filter_dbsnp:
        germline = []

    databse = get_database(species)
    exonic = get_exonic(exonic_file, species)
    translator = exonic["gene"].to_dict()

    # vcf = VCF(vcf_file, samples=list(info.index))
    vcf = VCF(vcf_file)
    cells = vcf.samples
    cells = pd.Series(cells, index=range(len(cells)), dtype=str, name="cells")
    n_total = 0
    n_multiallelic = 0
    n_germline = 0
    n_exonic = 0
    filtered = []
    for v in vcf():
        n_total += 1
        index = ""
        if len(v.ALT) != 1:
            n_multiallelic += 1
            continue
        if v.is_snp:
            index = f"{v.CHROM}.{v.start+1}.{v.REF}.{v.ALT[0]}"
        if not only_snp:
            if v.is_indel:
                a = v.REF
                b = v.ALT[0]
                if v.var_subtype == "del":
                    a = a[len(b) :]
                    b = "-"
                    index = f"{v.CHROM}.{v.start+2}.{a}.{b}"
                elif v.var_subtype == "ins":
                    a = "-"
                    b = b[len(a) :]
                    index = f"{v.CHROM}.{v.start+1}.{a}.{b}"
                if filter_indels_greater_than_2:
                    if len(a) > 2:
                        continue
                    if len(b) > 2:
                        continue
        if index in germline:
            n_germline += 1
            continue
        if only_exonic:
            if index not in translator:
                n_exonic += 1
                continue
        if filter_private:
            if sum(v.gt_types == 1) + sum(v.gt_types == 3) < 2:
                continue
        filtered.append(v)

    tsc.logg.info(f"{n_total} variants reported by GATK")
    tsc.logg.info(f"filtered {n_multiallelic} due to being multi allelic")
    tsc.logg.info(f"filtered {n_germline} due to being germline")
    tsc.logg.info(f"filtered {n_exonic} due to not being exonic")

    muts = []
    for v in filtered:
        if v.is_snp:
            index = f"{v.CHROM}.{v.start+1}.{v.REF}.{v.ALT[0]}"
        if v.is_indel:
            a = v.REF
            b = v.ALT[0]
            if v.var_subtype == "del":
                a = a[len(b) :]
                b = "-"
                index = f"{v.CHROM}.{v.start+2}.{a}.{b}"
            elif v.var_subtype == "ins":
                a = "-"
                b = b[len(a) :]
                index = f"{v.CHROM}.{v.start+1}.{a}.{b}"
        if only_exonic:
            mut = f"{translator[index]}.{index}"
        else:
            mut = index
        muts.append(mut)
    muts = pd.Series(muts, index=range(len(muts)), dtype=str, name="muts")

    if only_exonic:
        exonic = exonic.reset_index()
        exonic["muts"] = exonic["gene"] + "." + exonic["muts"]
        exonic = exonic[["muts", "kind", "amino_acid_change"]]
        exonic[["ensemble", "gene", "chrom", "position", "reference", "alteration"]] = (
            exonic["muts"].apply(tsc.ul.split_mut).apply(pd.Series)
        )
        exonic = exonic.set_index("muts")

    X = 3 * np.ones((len(cells), len(muts)), dtype=int)
    T = np.zeros((len(cells), len(muts)), dtype=int)
    V = np.zeros((len(cells), len(muts)), dtype=int)
    G = np.zeros((len(cells), len(muts)), dtype=int)

    for j, v in enumerate(filtered):
        T[:, j] = v.gt_ref_depths + v.gt_alt_depths
        V[:, j] = v.gt_alt_depths
        G[:, j] = v.gt_types

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(cells).set_index("cells"),
        var=exonic.loc[muts] if only_exonic else pd.DataFrame(index=muts),
        layers={"total": T, "mutant": V, "genotype": G},
    )
    if os.path.exists(f"{folderpath}/readstat.tsv"):
        df = pd.read_table(f"{folderpath}/readstat.tsv", index_col=0)
        adata.obs = pd.merge(
            adata.obs, df, how="left", right_index=True, left_index=True
        )

    tsc.logg.info(adata)
    return adata


def read_rsem(folderpath):
    """Read RSEM expression output to AnnData matrix.

    Parameters
    ----------
    folderpath : :obj:`str`
        Path to the folder contains the following files:

        1) `expr.count.tsv`
        2) `expr.fpkm.tsv`
        3) `expr.tpm.tsv`
        4) `readstat.tsv`

    Returns
    -------
    :class:`anndata.AnnData`
        The AnnData object which includes layers of FPKM and TPM values
        and count in .X values.
    """

    count = pd.read_table(f"{folderpath}/expr.count.tsv", index_col=0)
    tpm = pd.read_table(f"{folderpath}/expr.tpm.tsv", index_col=0)
    fpkm = pd.read_table(f"{folderpath}/expr.fpkm.tsv", index_col=0)
    info = pd.read_table(f"{folderpath}/readstat.tsv", index_col=0)
    info.index.name = "cells"

    adata = ad.AnnData(
        X=count.loc[info.index].values,
        obs=info,
        var=pd.DataFrame.from_dict({"genes": count.columns}).set_index("genes"),
        layers={"tpm": tpm.loc[info.index], "fpkm": fpkm.loc[info.index]},
    )

    tsc.logg.info(adata)
    return adata


def read_cnvkit(folderpath):
    data = {}
    for file in glob.glob(f"{folderpath}/*.cnr"):
        name = file.split("/")[-1].replace(".cnr", "")
        df = pd.read_table(file)
        df = df[(df.gene != "Antitarget")]
        df = df.dropna()
        data[name] = df

    n = len(data.keys())
    m = data[list(data.keys())[0]].shape[0]
    depth = np.zeros((n, m))
    log2 = np.zeros((n, m))
    weight = np.zeros((n, m))
    X = np.zeros((n, m))

    for i, cell in enumerate(data.keys()):
        depth[i, :] = data[cell].depth.values
        log2[i, :] = data[cell].log2.values
        weight[i, :] = data[cell].weight.values

    depth[-1, :] = np.zeros(m)
    log2[-1, :] = np.zeros(m)
    weight[-1, :] = np.zeros(m)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=list(data.keys())),
        var=data[list(data.keys())[0]][["chromosome", "start", "end", "gene"]],
        layers={"depth": depth, "log2": log2, "weight": weight},
    )

    cnv = np.zeros(adata.shape, dtype=int)
    tmp = adata.layers["log2"] < -1.1
    cnv[tmp == True] = 0
    tmp = (-1.1 <= adata.layers["log2"]) & (adata.layers["log2"] < -0.25)
    cnv[tmp == True] = 1
    tmp = (-0.25 <= adata.layers["log2"]) & (adata.layers["log2"] < 0.2)
    cnv[tmp == True] = 2
    tmp = (0.2 <= adata.layers["log2"]) & (adata.layers["log2"] < 0.7)
    cnv[tmp == True] = 3
    tmp = 0.7 <= adata.layers["log2"]
    cnv[tmp == True] = 4
    adata.layers["cnv"] = cnv

    abr = np.zeros(adata.shape, dtype=int)
    tmp = cnv < 2
    abr[tmp == True] = -1
    tmp = cnv == 2
    abr[tmp == True] = 0
    tmp = 2 < cnv
    abr[tmp == True] = 1
    adata.layers["abr"] = abr

    return adata


def read_bamreadcount(filepath):
    def get_nucleotides(line):
        chrom = line.split()[0]
        pos = line.split()[1]
        ref = line.split()[2].upper()
        depth = int(line.split()[3])
        char = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        char[ref] = int(line.split()[4].split(":")[1])
        char[line.split()[5].split(":")[0]] = int(line.split()[5].split(":")[1])
        char[line.split()[6].split(":")[0]] = int(line.split()[6].split(":")[1])
        char[line.split()[7].split(":")[0]] = int(line.split()[7].split(":")[1])
        char[line.split()[8].split(":")[0]] = int(line.split()[8].split(":")[1])
        nucleotides = " ".join([str(char[a]) for a in ["A", "C", "G", "T", "N"]])
        cov = sum([char[a] for a in ["A", "C", "G", "T", "N"]])
        sorted_char = sorted(char.items(), key=operator.itemgetter(1), reverse=True)
        return ref, char, cov

    data = []
    for file in glob.glob(f"{filepath}/*.txt"):
        with open(file) as fin:
            name = file.split("/")[-1].replace(".txt", "")
            for line in fin:
                line = line.strip()
                a = ".".join(line.split("\t")[:3])
                ref, char, cov = get_nucleotides(line)
                res = {"name": name, "mut": a, "ref": ref, "cov": cov}
                for k, v in char.items():
                    res[k] = v
                data.append(res)
    df = pd.DataFrame(data)
    df["index"] = df.apply(lambda x: f"{x['name']}.{x['mut']}", axis=1)

    return df


def read_sc_bulk_simulation(sc_file, bulk_file):
    """Read single-cell and bulk daat from Salem's simulation
    and merge them as an AnnData object.

    Parameters
    ----------
    sc_file : :obj:`str`
        Path to single-cell file.
    bulk_file : :obj:`str`
        Path to bulk file.

    Returns
    -------
    :class:`anndata.AnnData`
        An object contains VAF information in `.var`.
    """

    sc = tsc.io.read(sc_file)
    bulk = pd.read_table(bulk_file)
    bulk = bulk.set_index("ID")
    bulk[["trueNode", "trueVAF", "ID", "sampleIDs"]] = (
        bulk.INFO.str.split(";")
        .apply(pd.Series)[[0, 1, 2, 3]]
        .rename(columns={0: "trueNode", 1: "trueVAF", 2: "ID", 3: "sampleIDs"})
    )
    bulk = bulk.drop(["INFO"], axis=1)
    bulk.trueNode = bulk.trueNode.str.replace("trueNode=", "")
    bulk.trueVAF = bulk.trueVAF.str.replace("trueVAF=", "")
    bulk.ID = bulk.ID.str.replace("ID=", "")
    bulk.sampleIDs = bulk.sampleIDs.str.replace("sampleIDs=", "")
    adata = ad.AnnData(sc, dtype=int)
    adata.var = bulk

    return adata


def read_snpeff(filepath):
    """Read the VCF file annotated by SnpEff in multi-sample format.

    SnpEff was introduced in :cite:`SnpEff`.

    Parameters
    ----------
    filepath : :obj:`str`
        The path to the VCF file.

    Returns
    -------
    :class:`anndata.AnnData`
        The AnnData object which includes layers of mutant, total, genotype anc cnv.
    """

    vcf = VCF(filepath)
    info = vcf.get_header_type("ANN")["Description"].split(" | ")
    info[0] = "Allele"
    info[-1] = "ERRORS / WARNINGS / INFO"
    info = ["CHROM", "POS", "REF", "ALT", "START", "END"] + info

    cells = vcf.samples
    muts = []
    m_gen = []
    m_ref = []
    m_alt = []
    m_cna = []
    for var in VCF(filepath):
        if var.is_snp or var.is_indel:
            row = [var.CHROM, var.POS, var.REF, var.ALT, var.start + 1, var.end]
            ann = var.INFO.get("ANN").split(",")[0].split("|")
            row += ann
            m_gen.append(var.gt_types.tolist())
            m_ref.append(var.gt_ref_depths.tolist())
            m_alt.append(var.gt_alt_depths.tolist())
            muts.append(row)
            cnas_per_var = []
            for cell in cells:
                cna_tmp = var.INFO.get(f"cn|{cell}")
                if cna_tmp is not None:
                    cnas_per_var.append(int(cna_tmp))
                else:
                    cnas_per_var.append(-1)
            m_cna.append(cnas_per_var)

    m_gen = np.array(m_gen)
    m_ref = np.array(m_ref)
    m_alt = np.array(m_alt)
    m_cna = np.array(m_cna)
    muts = pd.DataFrame(muts, columns=info)
    muts.index = muts.index.map(lambda x: f"mut{x}")
    cells = pd.DataFrame(index=cells)

    adata = ad.AnnData(np.zeros((len(muts), len(cells))))
    adata.obs = muts
    adata.var = cells
    adata.layers["genotype"] = m_gen
    adata.layers["total"] = m_ref + m_alt
    adata.layers["mutant"] = m_alt
    adata.layers["cna"] = m_cna
    adata = adata.T
    return adata


def read_vep(filepath):
    """Read the VCF file annotated by VEP in multi-sample format.

    VEP was introduced in :cite:`VEP`.

    Parameters
    ----------
    filepath : :obj:`str`
        The path to the VCF file.

    Returns
    -------
    :class:`anndata.AnnData`
        The AnnData object which includes layers of mutant, total, genotype anc cnv.
    """

    vcf = VCF(filepath)
    info = vcf.get_header_type("CSQ")["Description"].split("|")
    info[0] = "Allele"
    info = ["CHROM", "POS", "REF", "ALT", "START", "END"] + info

    cells = vcf.samples
    muts = []
    m_gen = []
    m_ref = []
    m_alt = []
    m_cna = []
    for var in VCF(filepath):
        if var.is_snp or var.is_indel:
            row = [var.CHROM, var.POS, var.REF, var.ALT, var.start + 1, var.end]
            ann = var.INFO.get("CSQ").split(",")[0].split("|")
            row += ann
            m_gen.append(var.gt_types.tolist())
            m_ref.append(var.gt_ref_depths.tolist())
            m_alt.append(var.gt_alt_depths.tolist())
            muts.append(row)
            cnas_per_var = []
            for cell in cells:
                cna_tmp = var.INFO.get(f"cn|{cell}")
                if cna_tmp is not None:
                    cnas_per_var.append(int(cna_tmp))
                else:
                    cnas_per_var.append(-1)
            m_cna.append(cnas_per_var)

    m_gen = np.array(m_gen)
    m_ref = np.array(m_ref)
    m_alt = np.array(m_alt)
    m_cna = np.array(m_cna)
    muts = pd.DataFrame(muts, columns=info)
    muts.index = muts.index.map(lambda x: f"mut{x}")
    cells = pd.DataFrame(index=cells)

    adata = ad.AnnData(np.zeros((len(muts), len(cells))))
    adata.obs = muts
    adata.var = cells
    adata.layers["genotype"] = m_gen
    adata.layers["total"] = m_ref + m_alt
    adata.layers["mutant"] = m_alt
    adata.layers["cna"] = m_cna
    adata = adata.T
    return adata
