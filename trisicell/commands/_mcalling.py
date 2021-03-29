import csv
import glob
import gzip
import math
import os
import shutil
import subprocess

import click
import pandas as pd
import yaml
from natsort import natsorted

import trisicell as tsc

from .mcalling.s01indexing import run01
from .mcalling.s02mapping import run02
from .mcalling.s03indexing import run03
from .mcalling.s04mapping import run04
from .mcalling.s05calling import run05
from .mcalling.s06jointcalling import run06
from .mcalling.s07merging import run07
from .mcalling.s08annotating import run08
from .mcalling.s09expressing import run09
from .mcalling.s10velocitying import run10
from .mcalling.z01status import after01
from .mcalling.z02matrixing import after02
from .mcalling.z03clean import after03


def find_samples_readlength(config):
    def remove_pre_post(infile):
        if config["pairedend"]:
            return (
                os.path.basename(infile)
                .replace(config["infqpost1"], "")
                .replace(config["infqpost2"], "")
            )
        else:
            return os.path.basename(infile).replace(config["infqpost"], "")

    def get_readlength(filename):
        readlength = None
        with gzip.open(filename, "rb") as fin:
            i = 0
            for line in fin:
                readlength = len(line.strip()) - 1
                if i == 1:
                    break
                i += 1
        return readlength

    if os.path.exists(config["infq"]):
        infiles = glob.glob(f"{config['infq']}/*.fastq.gz")
        samples = natsorted(list(set(map(remove_pre_post, infiles))))
        readlength = get_readlength(infiles[0])
    else:
        infiles = glob.glob(f"{config['outdir']}/*/")
        samples = []
        for file in infiles:
            sample = file.split("/")[-2]
            if sample[0] != "_":
                samples.append(sample)
        samples = natsorted(samples)
        readlength = 73  # TODO: extract automatically
    return samples, readlength


@click.command(short_help="Mutation calling.")
@click.argument(
    "config_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.option(
    "--test",
    "-t",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Is in test mode.",
)
@click.option(
    "--status",
    "-s",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Check the status of runs.",
)
@click.option(
    "--build",
    "-b",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Build the final matrices of genotype and expression.",
)
@click.option(
    "--clean",
    "-c",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Cleaning the directory including removing intermediate files.",
)
def mcalling(config_file, test, status, build, clean):
    """Mutation calling pipeline.

    trisicell mcalling config.yml
    """

    with open(config_file) as fin:
        config = yaml.load(fin, Loader=yaml.FullLoader)

        samples, readlength = find_samples_readlength(config)
        config["samples"] = samples
        config["readlength"] = readlength

        tsc.logg.info(f"isRNA = {config['isrna']}")
        tsc.logg.info(config["samples"][:2], len(config["samples"]))
        tsc.logg.info(config["readlength"])

        tsc.ul.mkdir(config["tmpdir"] + "/cmd")
        tsc.ul.mkdir(config["tmpdir"] + "/log")
        tsc.ul.mkdir(config["outdir"])

        if config["buildver"] in ["hg19", "hg38"]:
            config["chroms"] = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        elif config["buildver"] == "mm10":
            config["chroms"] = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]

        if "infqpost1" in config and "infqpost2" in config:
            config["infqpost3"] = f"{config['infqpost1'].split('.')[0]}_val_1.fq.gz"
            config["infqpost4"] = f"{config['infqpost2'].split('.')[0]}_val_2.fq.gz"

        if "infqpost" in config:
            config["infqpost5"] = f"{config['infqpost']}_trimmed.fq.gz"

        #### STEPS
        if (not status) and (not build) and (not clean):
            if config["isrna"] == True:
                out01 = None
                if config["steps"]["s01indexing"]:
                    cmd01 = run01(config)
                    if not test:
                        out01 = subprocess.getoutput(cmd01)

                out02 = None
                if config["steps"]["s02mapping"]:
                    cmd02 = run02(config, out01)
                    if not test:
                        out02 = subprocess.getoutput(cmd02)

                out03 = None
                if config["steps"]["s03indexing"]:
                    cmd03 = run03(config, out02)
                    if not test:
                        out03 = subprocess.getoutput(cmd03)

                out04 = None
                if config["steps"]["s04mapping"]:
                    cmd04 = run04(config, out03)
                    if not test:
                        out04 = subprocess.getoutput(cmd04)

                out05 = None
                if config["steps"]["s05calling"]:
                    cmd05 = run05(config, out04)
                    if not test:
                        out05 = subprocess.getoutput(cmd05)

                out06 = None
                if config["steps"]["s06jointcalling"]:
                    cmd06 = run06(config, out05)
                    if not test:
                        out06 = subprocess.getoutput(cmd06)

                out07 = None
                if config["steps"]["s07merging"]:
                    cmd07 = run07(config, out06)
                    if not test:
                        out07 = subprocess.getoutput(cmd07)

                out08 = None
                if config["steps"]["s08annotating"]:
                    cmd08 = run08(config, out07)
                    if not test:
                        out08 = subprocess.getoutput(cmd08)

                out09 = None
                if config["steps"]["s09expressing"]:
                    cmd09 = run09(config, out05)
                    if not test:
                        out09 = subprocess.getoutput(cmd09)

                out10 = None
                if config["steps"]["s10velocitying"]:
                    cmd10 = run10(config, out04)
                    if not test:
                        out10 = subprocess.getoutput(cmd10)
            else:
                out02 = None
                if config["steps"]["s02mapping"]:
                    cmd02 = run02(config, None)
                    if not test:
                        out02 = subprocess.getoutput(cmd02)

                out04 = None
                if config["steps"]["s04mapping"]:
                    cmd04 = run04(config, out02)
                    if not test:
                        out04 = subprocess.getoutput(cmd04)

                out05 = None
                if config["steps"]["s05calling"]:
                    cmd05 = run05(config, out04)
                    if not test:
                        out05 = subprocess.getoutput(cmd05)

                out06 = None
                if config["steps"]["s06jointcalling"]:
                    cmd06 = run06(config, out05)
                    if not test:
                        out06 = subprocess.getoutput(cmd06)

                out07 = None
                if config["steps"]["s07merging"]:
                    cmd07 = run07(config, out06)
                    if not test:
                        out07 = subprocess.getoutput(cmd07)

                out08 = None
                if config["steps"]["s08annotating"]:
                    cmd08 = run08(config, out07)
                    if not test:
                        out08 = subprocess.getoutput(cmd08)
        else:
            if status:
                after01(config)
            if build:
                after02(config)
            if clean:
                after03(config)

    return None
