#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# ======================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 12, 2020
# Description: creating an initial reference index for STAR
# ======================================================================================

import os

import pandas as pd

from trisicell.ul._servers import cmd, write_cmds_get_main


def run01(config):
    def cmds():
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/_indexing/1"])
        cmds += cmd(["module load STAR/2.7.3a"])
        cmds += cmd(
            [
                "STAR",
                "--runMode genomeGenerate",
                f"--genomeDir {config['outdir']}/_indexing/1",
                f"--genomeFastaFiles {config['ref']}",
                f"--sjdbGTFfile {config['annot']}",
                f"--sjdbOverhang {config['readlength']}",
                "--runThreadN 32",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["cmd"] = [cmds()]

    jobname = os.path.basename(__file__)[:-3]
    cmdmain = write_cmds_get_main(
        df,
        jobname,
        config["time01"],
        config["mem01"],
        "",
        32,
        config["email"],
        config["tmpdir"],
        None,
    )
    return cmdmain
