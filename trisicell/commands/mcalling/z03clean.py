#!/usr/bin/env python

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 14, 2020
# Description: cleaning
# =========================================================================================

import subprocess

from trisicell.ul._servers import *


def after03(config):
    # subprocess.getoutput(f'rm -rf {config['outdir']}/_annotating')
    subprocess.getoutput(f"rm -rf {config['outdir']}/_calling")
    subprocess.getoutput(f"rm -rf {config['outdir']}/_indexing")

    subprocess.getoutput(
        f"find {config['outdir']} -name \\*Aligned.out.sam -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*Aligned.sortedByCoord.out.bam -exec rm -rf"
        " {} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*Aligned.toTranscriptome.out.bam -exec rm -rf"
        " {} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*dedupped.bam -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*expr.genes.results -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*expr.isoforms.results -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*expr.stat -exec rm -rf {{}} \\;"
    )
    # subprocess.getoutput(f'find {config['outdir']} -name \\*HaplotypeCaller.g.vcf -exec rm -rf {{}} \\;')
    # subprocess.getoutput(f'find {config['outdir']} -name \\*HaplotypeCaller.g.vcf.idx -exec rm -rf {{}} \\;')
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*indel.intervals -exec rm -rf {{}} \\;"
    )
    # subprocess.getoutput(f'find {config['outdir']} -name \\*Log.final.out -exec rm -rf {{}} \\;')
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*Log.out -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*Log.progress.out -exec rm -rf {{}} \\;"
    )
    # subprocess.getoutput(f'find {config['outdir']} -name \\*output.bai -exec rm -rf {{}} \\;')
    # subprocess.getoutput(f'find {config['outdir']} -name \\*output.bam -exec rm -rf {{}} \\;')
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*output.metrics -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*ReadsPerGene.out.tab -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*realign.bai -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*realign.bam -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*recal.table -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*rg_added_sorted.bai -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*rg_added_sorted.bam -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*SJ.out.tab -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*split.bai -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*split.bam -exec rm -rf {{}} \\;"
    )
    subprocess.getoutput(
        f"find {config['outdir']} -name \\*_STARgenome -exec rm -rf {{}} \\;"
    )

    if config["dotrimming"]:
        subprocess.getoutput(
            f"find {config['outdir']} -name \\*_fastqc.html -exec rm -rf {{}} \\;"
        )
        subprocess.getoutput(
            f"find {config['outdir']} -name \\*_report.txt -exec rm -rf {{}} \\;"
        )
        subprocess.getoutput(
            f"find {config['outdir']} -name \\*_fastqc.zip -exec rm -rf {{}} \\;"
        )
        subprocess.getoutput(
            f"find {config['outdir']} -name \\*.fq.gz -exec rm -rf {{}} \\;"
        )
