#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 14, 2020
# Description: mapping for the first pass of the STAR
# =========================================================================================

from trisicell.ul._servers import *


def run02(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if config["dotrimming"] == True:
            if config["pairedend"] == True:
                cmds += cmd(
                    [
                        f"trim_galore",
                        f"--gzip",
                        f"--length 30",
                        f"--fastqc",
                        f"--output_dir {config['outdir']}/{sample}",
                        f"--paired",
                        f"{config['infq']}/{config['infqpre1']}{sample}{config['infqpost1']} {config['infq']}/{config['infqpre2']}{sample}{config['infqpost2']}",
                    ]
                )
                if config["isrna"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            f"--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['outdir']}/{sample}/{config['infqpre1']}{sample}{config['infqpost3']} {config['outdir']}/{sample}/{config['infqpre2']}{sample}{config['infqpost4']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"bwa mem",
                            f"-t 8",
                            f"{config['ref']}",
                            f"{config['outdir']}/{sample}/{config['infqpre1']}{sample}{config['infqpost3']} {config['outdir']}/{sample}/{config['infqpre2']}{sample}{config['infqpost4']}",
                            f"|",
                            f"samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"samtools index",
                            f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
            else:
                cmds += cmd(
                    [
                        f"trim_galore",
                        f"--gzip",
                        f"--length 30",
                        f"--fastqc",
                        f"--output_dir {config['outdir']}/{sample}",
                        f"{config['infq']}/{config['infqpre']}{sample}{config['infqpost']}",
                    ]
                )
                if config["isrna"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            f"--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['outdir']}/{sample}/{config['infqpre']}{sample}{config['infqpost5']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"bwa mem",
                            f"-t 8",
                            f"{config['ref']}",
                            f"{config['outdir']}/{sample}/{config['infqpre']}{sample}{config['infqpost5']}",
                            f"|",
                            f"samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"samtools index",
                            f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
        else:
            if config["pairedend"] == True:
                if config["isrna"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            f"--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['infq']}/{config['infqpre1']}{sample}{config['infqpost1']} {config['infq']}/{config['infqpre2']}{sample}{config['infqpost2']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"bwa mem",
                            f"-t 8",
                            f"{config['ref']}",
                            f"{config['infq']}/{config['infqpre1']}{sample}{config['infqpost1']} {config['infq']}/{config['infqpre2']}{sample}{config['infqpost2']}",
                            f"|",
                            f"samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"samtools index",
                            f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
            else:
                if config["isrna"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/1",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                            f"--readFilesCommand zcat",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--readFilesIn {config['infq']}/{config['infqpre']}{sample}{config['infqpost']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/Aligned.out.sam",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"bwa mem",
                            f"-t 8",
                            f"{config['ref']}",
                            f"{config['infq']}/{config['infqpre']}{sample}{config['infqpost']}",
                            f"|",
                            f"samtools sort -@8 -m 10000000000",
                            f"-o {config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam -",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"samtools index",
                            f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                        ]
                    )
        cmds += cmd([f"echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    df["sample"] = config["samples"]
    df["cmd"] = df.apply(lambda x: cmds(x["sample"]), axis=1)

    jobname = os.path.basename(__file__)[:-3]
    if config["isrna"] == True:
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
            "bwa/0.7.17,samtools/1.11,trimmomatic/0.39,trimgalore/0.6.5,cutadapt/2.9,fastqc/0.11.9",
            8,
            config["email"],
            config["tmpdir"],
            afterok,
        )
    return cmdmain
