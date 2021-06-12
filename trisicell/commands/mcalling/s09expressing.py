#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Jan 1, 2021
# Description: get the expression using RSEM
# =========================================================================================

from trisicell.ul._servers import *


def run09(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if config["pairedend"]:
            cmds += cmd(
                [
                    f"rsem-calculate-expression",
                    f"--bam",
                    f"--no-bam-output",
                    f"--paired-end",
                    f"--estimate-rspd",
                    f"--append-names",
                    f"{config['outdir']}/{sample}/Aligned.toTranscriptome.out.bam",
                    f"{config['rsemdb']}",
                    f"{config['outdir']}/{sample}/expr",
                ]
            )
        else:
            cmds += cmd(
                [
                    f"rsem-calculate-expression",
                    f"--bam",
                    f"--no-bam-output",
                    f"--estimate-rspd",
                    f"--append-names",
                    f"{config['outdir']}/{sample}/Aligned.toTranscriptome.out.bam",
                    f"{config['rsemdb']}",
                    f"{config['outdir']}/{sample}/expr",
                ]
            )
        cmds += cmd(
            [
                f"rm -rf",
                f"{config['outdir']}/{sample}/expr.stat",
                f"{config['outdir']}/{sample}/expr.isoforms.results",
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
        config["time09"],
        config["mem09"],
        "STAR/2.7.3a,rsem/1.3.2",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
