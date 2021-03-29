#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: May 17, 2020
# Description: calling mutations as a joint cohort by GATK GenotypeGVCFs for each chr
# =========================================================================================

from trisicell.ul._servers import *


def run06(config, afterok):
    def cmds(chrom):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/_calling"])
        if len(config["normals"]) > 0:
            files = " -V ".join(
                [
                    f"{config['outdir']}/{s}/MuTect2.g.vcf"
                    for s in config["samples"]
                    if s not in config["normals"]
                ]
            )
        else:
            files = " -V ".join(
                [
                    f"{config['outdir']}/{s}/HaplotypeCaller.g.vcf"
                    for s in config["samples"]
                ]
            )
        cmds += cmd(
            [
                f"{config['java06']} {config['GATK']}",
                f"-T GenotypeGVCFs",
                f"-R {config['ref']}",
                f"-V {files}",
                f"-o {config['outdir']}/_calling/jointcalls.{chrom}.g.vcf",
                f"-nt {config['threads06']}",
                f"-L {chrom}",
                f"--disable_auto_index_creation_and_locking_when_reading_rods",
            ]
        )
        cmds += cmd([f"echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["chrom"] = config["chroms"]
    df["cmd"] = df.apply(lambda x: cmds(x["chrom"]), axis=1)

    jobname = os.path.basename(__file__)[:-3]
    cmdmain = write_cmds_get_main(
        df,
        jobname,
        config["time06"],
        config["mem06"],
        "python/3.7",
        config["threads06"],
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
