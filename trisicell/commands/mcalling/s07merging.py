#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 03, 2020
# Description: merging joint mutations called
# =========================================================================================

from trisicell.ul._servers import *


def run07(config, afterok):
    def cmds():
        cmds = ""
        cmds += cmd([f"module load bcftools/1.9"])
        files = " ".join(
            [
                f"{config['outdir']}/_calling/jointcalls.{chrom}.g.vcf"
                for chrom in config["chroms"]
            ]
        )
        cmds += cmd(
            [
                f"bcftools concat",
                f"-o {config['outdir']}/_calling/jointcalls.g.vcf",
                f"{files}",
            ]
        )
        cmds += cmd([f"echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["cmd"] = [cmds()]

    jobname = os.path.basename(__file__)[:-3]
    cmdmain = write_cmds_get_main(
        df,
        jobname,
        config["time07"],
        config["mem07"],
        "",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
