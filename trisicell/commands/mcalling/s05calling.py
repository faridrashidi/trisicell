#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: May 16, 2020
# Description: calling mutations by GATK HaplotypeCaller
# =========================================================================================

from trisicell.ul._servers import *


def run05(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if len(config["normals"]) > 0:
            if sample not in config["normals"]:
                cmds += cmd(
                    [
                        f"{config['java05']} {config['GATK']}",
                        f"-T MuTect2",
                        f"-R {config['ref']}",
                        f"-I:tumor {config['outdir']}/{sample}/output.bam",
                        "-I:normal"
                        f" {config['outdir']}/{config['normals'][0]}/output.bam",
                        f"-o {config['outdir']}/{sample}/MuTect2.g.vcf",
                        f"-dontUseSoftClippedBases",
                        f"-stand_call_conf 20",
                        f"--dbsnp {config['hgDBSNPs138']}"
                        if buildver == "hg19"
                        else f"--dbsnp {config['mmDBSNPs142']}",
                        f"-ERC GVCF",
                    ]
                )
        else:
            cmds += cmd(
                [
                    f"{config['java05']} {config['GATK']}",
                    f"-T HaplotypeCaller",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/output.bam",
                    f"-o {config['outdir']}/{sample}/HaplotypeCaller.g.vcf",
                    f"-dontUseSoftClippedBases",
                    f"-stand_call_conf 20",
                    f"--dbsnp {config['hgDBSNPs138']}"
                    if config["buildver"] == "hg19"
                    else f"--dbsnp {config['mmDBSNPs142']}",
                    f"-ERC GVCF",
                ]
            )
        cmds += cmd([f"echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["sample"] = config["samples"]
    df["cmd"] = df.apply(lambda x: cmds(x["sample"]), axis=1)

    jobname = os.path.basename(__file__)[:-3]
    cmdmain = write_cmds_get_main(
        df,
        jobname,
        config["time05"],
        config["mem05"],
        "python/3.7",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
