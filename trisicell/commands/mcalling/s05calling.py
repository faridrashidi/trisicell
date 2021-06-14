#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# ======================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: May 16, 2020
# Description: calling mutations by GATK HaplotypeCaller
# ======================================================================================

import os

import pandas as pd

from trisicell.ul._servers import cmd, write_cmds_get_main


def run05(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if len(config["normals"]) > 0:
            if sample not in config["normals"]:
                cmds += cmd(
                    [
                        f"{config['java05']} {config['GATK']}",
                        "-T MuTect2",
                        f"-R {config['ref']}",
                        f"-I:tumor {config['outdir']}/{sample}/output.bam",
                        "-I:normal"
                        f" {config['outdir']}/{config['normals'][0]}/output.bam",
                        f"-o {config['outdir']}/{sample}/MuTect2.g.vcf",
                        "-dontUseSoftClippedBases",
                        "-stand_call_conf 20",
                        f"--dbsnp {config['hgDBSNPs138']}"
                        if config["buildver"] == "hg19"
                        else f"--dbsnp {config['mmDBSNPs142']}",
                        "-ERC GVCF",
                    ]
                )
        else:
            cmds += cmd(
                [
                    f"{config['java05']} {config['GATK']}",
                    "-T HaplotypeCaller",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/output.bam",
                    f"-o {config['outdir']}/{sample}/HaplotypeCaller.g.vcf",
                    "-dontUseSoftClippedBases",
                    "-stand_call_conf 20",
                    f"--dbsnp {config['hgDBSNPs138']}"
                    if config["buildver"] == "hg19"
                    else f"--dbsnp {config['mmDBSNPs142']}",
                    "-ERC GVCF",
                ]
            )
        cmds += cmd(["echo Done!"], islast=True)
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
