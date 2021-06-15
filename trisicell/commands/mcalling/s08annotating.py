#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# ======================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: May 16, 2020
# Description: annotating mutations by ANNOVAR to find the missense ones
# ======================================================================================

import os

import pandas as pd

from trisicell.ul._servers import cmd, write_cmds_get_main


def run08(config, afterok):
    def cmds():
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/_annotating"])
        cmds += cmd([f"mkdir -p {config['outdir']}/_annotating/annovar_pre"])
        cmds += cmd([f"mkdir -p {config['outdir']}/_annotating/annovar_post"])
        cmds += cmd([f"mkdir -p {config['outdir']}/_annotating/annovar_final"])
        cmds += cmd(
            [
                f"{config['ANNOVAR']}/convert2annovar.pl",
                "--format vcf4",
                f"{config['outdir']}/_calling/jointcalls.g.vcf",
                f"--outfile {config['outdir']}/_annotating/annovar_pre/pre",
                "--includeinfo",
                "--allsample",
            ]
        )
        for file in [
            f"{config['outdir']}/_annotating/annovar_pre/pre.{s}.avinput"
            for s in config["samples"]
        ]:
            cmds += cmd(
                [
                    f"{config['ANNOVAR']}/annotate_variation.pl",
                    f"--outfile {config['outdir']}/_annotating/annovar_post/"
                    + f"post.{os.path.basename(file)[4:]}",
                    f"--dbtype {config['dbtype']}",
                    f"--buildver {config['buildver']}",
                    f"{file}",
                    f"{config['ANNOVAR']}/{config['whichdb']}/",
                ]
            )
        for file in [
            f"{config['outdir']}/_annotating/annovar_post/"
            + f"post.{s}.avinput.exonic_variant_function"
            for s in config["samples"]
        ]:
            cmds += cmd(
                [
                    'awk -F \'\\t\' \'$2!="unknown" {split($3,a,":");',
                    'print a[2]"_"a[1]"\\t"$4"\\t"$5"\\t"$7"\\t"$8"\\t"$2"\\t"$3}\'',
                    file,
                    f"> {config['outdir']}/_annotating/annovar_final/"
                    + f"{os.path.basename(file)[5:-32]}.txt",
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
        config["time08"],
        config["mem08"],
        "",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
