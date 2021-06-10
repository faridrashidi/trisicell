#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Dec 31, 2020
# Description: using velocyto for RNA velocity
# =========================================================================================

from trisicell.ul._servers import *


def run10(config, afterok):
    def cmds():
        cmds = ""
        cmds += cmd([f"module load velocyto/0.17"])
        cmds += cmd(
            [
                f"velocyto run-smartseq2",
                f"-o {config['outdir']}/_velocyto",
                f"-m {config['velo']}",
                f"-e velocyto",
                f"{config['outdir']}/*/*.align.bam",
                f"{config['annot']}",
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
        config["time10"],
        config["mem10"],
        "",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
