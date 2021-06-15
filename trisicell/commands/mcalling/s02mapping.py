#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# ======================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 14, 2020
# Description: mapping for the first pass of the STAR
# ======================================================================================

import os

import pandas as pd

from trisicell.ul._servers import cmd, write_cmds_get_main


def run02(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if config["dotrimming"]:
            if config["pairedend"]:
                cmds += cmd(
                    [
                        "trim_galore",
                        "--gzip",
                        "--length 30",
                        "--fastqc",
                        f"--output_dir {config['outdir']}/{sample}",
                        "--paired",
                        f"{config['infq']}/{config['infqpre1']}"
                        + f"{sample}{config['infqpost1']}"
                        f" {config['infq']}/{config['infqpre2']}"
                        + f"{sample}{config['infqpost2']}",
                    ]
                )
                if config["isrna"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            "--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost3']}"
                            f" {config['outdir']}/{sample}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost4']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "bwa mem",
                            "-t 8",
                            f"{config['ref']}",
                            f"{config['outdir']}/{sample}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost3']}"
                            f" {config['outdir']}/{sample}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost4']}",
                            "|",
                            "samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            "samtools index",
                            f"{config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
            else:
                cmds += cmd(
                    [
                        "trim_galore",
                        "--gzip",
                        "--length 30",
                        "--fastqc",
                        f"--output_dir {config['outdir']}/{sample}",
                        f"{config['infq']}/{config['infqpre']}"
                        + f"{sample}{config['infqpost']}",
                    ]
                )
                if config["isrna"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            "--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre']}"
                            + f"{sample}{config['infqpost5']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "bwa mem",
                            "-t 8",
                            f"{config['ref']}",
                            f"{config['outdir']}/{sample}/{config['infqpre']}"
                            + f"{sample}{config['infqpost5']}",
                            "|",
                            "samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            "samtools index",
                            f"{config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
        else:
            if config["pairedend"]:
                if config["isrna"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            "--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['infq']}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost1']}"
                            f" {config['infq']}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost2']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "bwa mem",
                            "-t 8",
                            f"{config['ref']}",
                            f"{config['infq']}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost1']}"
                            f" {config['infq']}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost2']}",
                            "|",
                            "samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            "samtools index",
                            f"{config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
            else:
                if config["isrna"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            "--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['infq']}/{config['infqpre']}"
                            + f"{sample}{config['infqpost']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "bwa mem",
                            "-t 8",
                            f"{config['ref']}",
                            f"{config['infq']}/{config['infqpre']}"
                            + f"{sample}{config['infqpost']}",
                            "|",
                            "samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            "samtools index",
                            f"{config['outdir']}/"
                            + f"{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["sample"] = config["samples"]
    df["cmd"] = df.apply(lambda x: cmds(x["sample"]), axis=1)

    jobname = os.path.basename(__file__)[:-3]
    if config["isrna"]:
        cmdmain = write_cmds_get_main(
            df,
            jobname,
            config["time02"],
            config["mem02"],
            "STAR/2.7.3a,trimmomatic/0.39,trimgalore/0.6.5,cutadapt/2.9,fastqc/0.11.9",
            1,
            config["email"],
            config["tmpdir"],
            afterok,
        )
    else:
        cmdmain = write_cmds_get_main(
            df,
            jobname,
            config["time02"],
            config["mem02"],
            "bwa/0.7.17,samtools/1.11,trimmomatic/0.39,trimgalore/0.6.5,cutadapt/2.9,"
            + "fastqc/0.11.9",
            8,
            config["email"],
            config["tmpdir"],
            afterok,
        )
    return cmdmain
