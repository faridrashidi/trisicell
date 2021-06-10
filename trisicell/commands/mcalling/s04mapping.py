#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: July 11, 2020
# Description: mapping for the second pass of the STAR + GATK best practices
# =========================================================================================

from trisicell.ul._servers import *


def run04(config, afterok):
    def cmds(sample):
        cmds = ""
        cmds += cmd([f"mkdir -p {config['outdir']}/{sample}"])
        if config["isrna"] == True:
            if config["dotrimming"] == True:
                if config["pairedend"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            f"--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre1']}{sample}{config['infqpost3']}"
                            f" {config['outdir']}/{sample}/{config['infqpre2']}{sample}{config['infqpost4']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--limitBAMsortRAM 30000000000",
                            f"--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/{config['infqpre1']}{sample}{config['infqpost3']}",
                            f"{config['outdir']}/{sample}/{config['infqpre2']}{sample}{config['infqpost4']}",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            f"--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['outdir']}/{sample}/{config['infqpre']}{sample}{config['infqpost5']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--limitBAMsortRAM 30000000000",
                            f"--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                    cmds += cmd(
                        [
                            f"rm -rf",
                            f"{config['outdir']}/{sample}/{config['infqpre']}{sample}{config['infqpost5']}",
                        ]
                    )
            else:
                if config["pairedend"] == True:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            f"--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['infq']}/{config['infqpre1']}{sample}{config['infqpost1']}"
                            f" {config['infq']}/{config['infqpre2']}{sample}{config['infqpost2']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--limitBAMsortRAM 30000000000",
                            f"--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
                else:
                    cmds += cmd(
                        [
                            f"STAR",
                            f"--runMode alignReads",
                            f"--genomeDir {config['outdir']}/_indexing/2",
                            f"--readFilesCommand zcat",
                            "--readFilesIn"
                            f" {config['infq']}/{config['infqpre']}{sample}{config['infqpost']}",
                            f"--outFileNamePrefix {config['outdir']}/{sample}/",
                            f"--limitBAMsortRAM 30000000000",
                            f"--outSAMtype BAM SortedByCoordinate",
                            f"--sjdbGTFfile {config['annot']}",
                            f"--outFilterMultimapNmax 1",
                            f"--outSAMunmapped None",
                            f"--quantMode TranscriptomeSAM GeneCounts",
                            f"--runThreadN 1",
                            f"--sjdbOverhang {config['readlength']}",
                        ]
                    )
        cmds += cmd(
            [
                f"rm -rf",
                f"{config['outdir']}/{sample}/Aligned.out.sam",
            ]
        )

        cmds += cmd(
            [
                f"{config['java04']} {config['PICARD']}",
                f"AddOrReplaceReadGroups",
                f"INPUT={config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                f"OUTPUT={config['outdir']}/{sample}/dedupped.bam",
                f"SORT_ORDER=coordinate",
                f"RGID={sample}",
                f"RGLB=trancriptome" if config["isrna"] == True else f"RGLB=genome",
                f"RGPL=ILLUMINA",
                f"RGPU=machine",
                f"RGSM={sample}",
            ]
        )

        cmds += cmd(
            [
                f"{config['java04']} {config['PICARD']}",
                f"MarkDuplicates",
                f"INPUT={config['outdir']}/{sample}/dedupped.bam",
                f"OUTPUT={config['outdir']}/{sample}/rg_added_sorted.bam",
                f"METRICS_FILE={config['outdir']}/{sample}/output.metrics",
                f"VALIDATION_STRINGENCY=SILENT",
                f"CREATE_INDEX=true",
            ]
        )
        cmds += cmd(
            [
                f"rm -rf",
                f"{config['outdir']}/{sample}/dedupped.bam",
            ]
        )

        if config["isrna"] == True:
            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    f"-T SplitNCigarReads",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/split.bam",
                    f"-rf ReassignOneMappingQuality",
                    f"-RMQF 255",
                    f"-RMQT 60",
                    f"-U ALLOW_N_CIGAR_READS",
                ]
            )

        if config["buildver"] == "hg19":
            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    f"-T RealignerTargetCreator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"] == True
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['hgKGIndels']}",
                    f"-known {config['hgMillsIndels']}",
                    f"-U ALLOW_SEQ_DICT_INCOMPATIBILITY",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    f"-T IndelRealigner",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"] == True
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
                    f"-T BaseRecalibrator",
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
                    f"-T RealignerTargetCreator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"] == True
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['mmDBIndels142']}",
                    f"-U ALLOW_SEQ_DICT_INCOMPATIBILITY",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    f"-T IndelRealigner",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/split.bam"
                    if config["isrna"] == True
                    else f"-I {config['outdir']}/{sample}/rg_added_sorted.bam",
                    f"-o {config['outdir']}/{sample}/realign.bam",
                    f"-targetIntervals {config['outdir']}/{sample}/indel.intervals",
                    f"-known {config['mmDBIndels142']}",
                ]
            )

            cmds += cmd(
                [
                    f"{config['java04']} {config['GATK']}",
                    f"-T BaseRecalibrator",
                    f"-R {config['ref']}",
                    f"-I {config['outdir']}/{sample}/realign.bam",
                    f"-o {config['outdir']}/{sample}/recal.table",
                    f"-knownSites {config['mmDBSNPs142']}",
                    f"-knownSites {config['mmDBIndels142']}",
                ]
            )
        cmds += cmd(
            [
                f"rm -rf",
                f"{config['outdir']}/{sample}/rg_added_sorted.bam",
                f"{config['outdir']}/{sample}/rg_added_sorted.bai",
                f"{config['outdir']}/{sample}/split.bam",
                f"{config['outdir']}/{sample}/split.bai",
            ]
        )
        cmds += cmd(
            [
                f"{config['java04']} {config['GATK']}",
                f"-T PrintReads",
                f"-R {config['ref']}",
                f"-I {config['outdir']}/{sample}/realign.bam",
                f"-o {config['outdir']}/{sample}/output.bam",
                f"-BQSR {config['outdir']}/{sample}/recal.table",
            ]
        )
        cmds += cmd(
            [
                f"rm -rf",
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
                f"mv",
                f"{config['outdir']}/{sample}/Aligned.sortedByCoord.out.bam",
                f"{config['outdir']}/{sample}/{sample}.align.bam",
            ]
        )
        cmds += cmd(
            [
                f"samtools index",
                f"{config['outdir']}/{sample}/{sample}.align.bam",
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
        config["time04"],
        config["mem04"],
        "STAR/2.7.3a,samtools/1.11",
        1,
        config["email"],
        config["tmpdir"],
        afterok,
    )
    return cmdmain
