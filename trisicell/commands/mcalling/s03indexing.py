#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: May 16, 2020
# Description: creating a new reference index based on the first pass run of the STAR
# =========================================================================================

from trisicell.ul._servers import *


def run03(config, afterok):
    def cmds():
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/_indexing/2"])
        cmds += cmd([f"module load STAR/2.7.3a"])
        files = " ".join(
            [f"{config['outdir']}/{s}/SJ.out.tab" for s in config["samples"]]
        )
        cmds += cmd(
            [
                f"STAR",
                f"--runMode genomeGenerate",
                f"--genomeDir {config['outdir']}/_indexing/2",
                f"--genomeFastaFiles {config['ref']}",
                f"--sjdbGTFfile {config['annot']}",
                f"--sjdbFileChrStartEnd {files}",
                f"--sjdbOverhang {config['readlength']}",
                f"--runThreadN 32",
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
        config["time03"],
        config["mem03"],
        "",
        32,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
