#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# ======================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: July 11, 2020
# Description: mapping for the second pass of the STAR + GATK best practices
# ======================================================================================

import os

import pandas as pd

from trisicell.ul._servers import cmd, write_cmds_get_main


def run04(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if config["isrna"]:
            if config["dotrimming"]:
                if config["pairedend"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            "--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost3']}"
                            f" {config['outdir']}/{sample}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost4']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--limitBAMsortRAM 30000000000",
                            "--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost3']}",
                            f"{config['outdir']}/{sample}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost4']}",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            "--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre']}"
                            + f"{sample}{config['infqpost5']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--limitBAMsortRAM 30000000000",
                            "--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            "rm -rf",
                            f"{config['outdir']}/{sample}/{config['infqpre']}"
                            + f"{sample}{config['infqpost5']}",
                        ]
                    )
            else:
                if config["pairedend"]:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            "--readFilesCommand zcat",
                            f"--readFilesIn {config['infq']}/{config['infqpre1']}"
                            + f"{sample}{config['infqpost1']}"
                            f" {config['infq']}/{config['infqpre2']}"
                            + f"{sample}{config['infqpost2']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--limitBAMsortRAM 30000000000",
                            "--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            "STAR",
                            "--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            "--readFilesCommand zcat",
                            f"--readFilesIn {config['infq']}/{config['infqpre']}"
                            + f"{sample}{config['infqpost']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            "--limitBAMsortRAM 30000000000",
                            "--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            "--outFilterMultimapNmax 1",
                            "--outSAMunmapped None",
                            "--quantMode TranscriptomeSAM GeneCounts",
                            "--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
        cmds += cmd(
            [
                "rm -rf",
                f"{config['outdir']}/{sample}/Aligned.out.sam",
            ]
        )

        cmds += cmd(
            [
                f"{config['java04']} {config['PICARD']}",
                "AddOrReplaceReadGroups",
                f"INPUT={config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                f"OUTPUT={config['outdir']}/{sample}/dedupped.bam",
                "SORT_ORDER=coordinate",
                f"RGID={sample}",
                "RGLB=trancriptome" if config["isrna"] else "RGLB=genome",
                "RGPL=ILLUMINA",
                "RGPU=machine",
                f"RGSM={sample}",
            ]
        )

        cmds += cmd(
            [
                f"{config['java04']} {config['PICARD']}",
                "MarkDuplicates",
                f"INPUT={config['outdir']}/{sample}/dedupped.bam",
                f"OUTPUT={config['outdir']}/{sample}/rg_added_sorted.bam",
                f"METRICS_FILE={config['outdir']}/{sample}/output.metrics",
                "VALIDATION_STRINGENCY=SILENT",
                "CREATE_INDEX=true",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{config['outdir']}/{sample}/dedupped.bam",
            ]
        )

        if config["isrna"]:
            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T SplitNCigarReads",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/split.bam",
                    "-rf ReassignOneMappingQuality",
                    "-RMQF 255",
                    "-RMQT 60",
                    "-U ALLOW_N_CIGAR_READS",
                ]
            )

        if config["buildver"] == "hg19":
            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T RealignerTargetCreator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"]
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['hgKGIndels']}",
                    f"-known {config['hgMillsIndels']}",
                    "-U ALLOW_SEQ_DICT_INCOMPATIBILITY",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T IndelRealigner",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"]
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/realign.bam",
                    f"-targetIntervals {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['hgKGIndels']}",
                    f"-known {config['hgMillsIndels']}",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T BaseRecalibrator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/realign.bam",
                    f"-o {config['outdir']}/{sample}/recal.table",
                    f"-knownSites {config['hgDBSNPs138']}",
                    f"-knownSites {config['hgKGIndels']}",
                    f"-knownSites {config['hgMillsIndels']}",
                ]
            )
        else:
            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T RealignerTargetCreator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"]
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['mmDBIndels142']}",
                    "-U ALLOW_SEQ_DICT_INCOMPATIBILITY",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T IndelRealigner",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"]
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/realign.bam",
                    f"-targetIntervals {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['mmDBIndels142']}",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    "-T BaseRecalibrator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/realign.bam",
                    f"-o {config['outdir']}/{sample}/recal.table",
                    f"-knownSites {config['mmDBSNPs142']}",
                    f"-knownSites {config['mmDBIndels142']}",
                ]
            )
        cmds += cmd(
            [
                "rm -rf",
                f"{config['outdir']}/{sample}/rg_added_sorted.bam",
                f"{config['outdir']}/{sample}/rg_added_sorted.bai",
                f"{config['outdir']}/{sample}/split.bam",
                f"{config['outdir']}/{sample}/split.bai",
            ]
        )
        cmds += cmd(
            [
                f"{config['java04']} {config['GATK']}",
                "-T PrintReads",
                f"-R {config['ref']}",
                f"-I {config['outdir']}/{sample}/realign.bam",
                f"-o {config['outdir']}/{sample}/output.bam",
                f"-BQSR {config['outdir']}/{sample}/recal.table",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{config['outdir']}/{sample}/realign.bam",
                f"{config['outdir']}/{sample}/realign.bai",
                f"{config['outdir']}/{sample}/Log.progress.out",
                f"{config['outdir']}/{sample}/Log.out",
                f"{config['outdir']}/{sample}/indel.intervals",
                f"{config['outdir']}/{sample}/output.metrics",
                f"{config['outdir']}/{sample}/ReadsPerGene.out.tab",
                f"{config['outdir']}/{sample}/recal.table",
                f"{config['outdir']}/{sample}/SJ.out.tab",
                f"{config['outdir']}/{sample}/_STARgenome",
                # f'{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam',
            ]
        )
        cmds += cmd(
            [
                "mv",
                f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                f"{config['outdir']}/{sample}/{sample}.align.bam",
            ]
        )
        cmds += cmd(
            [
                "samtools index",
                f"{config['outdir']}/{sample}/{sample}.align.bam",
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
        config["time04"],
        config["mem04"],
        "STAR/2.7.3a,samtools/1.11",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
