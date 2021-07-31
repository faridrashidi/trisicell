import os

import click
import pandas as pd

import trisicell as tsc
from trisicell.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run SRA.")
@click.argument(
    "sra_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "out_dir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.option(
    "--check",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Check the unfinished jobs on biowulf.",
)
def sra(sra_file, out_dir, check):
    """Run SRA.

    trisicell sra path/to/SraRunTable.csv path/to/output/directory
    """

    def cmds(srr_id, name, layout):
        cmds = ""
        if layout == "PAIRED":
            cmds += cmd(
                [
                    "fastq-dump",
                    "--split-files",
                    f"--outdir {out_dir}",
                    f"{srr_id}",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{out_dir}/{srr_id}_1.fastq",
                    f"{out_dir}/{name}_1.fastq",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{out_dir}/{srr_id}_2.fastq",
                    f"{out_dir}/{name}_2.fastq",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{out_dir}/{name}_1.fastq",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{out_dir}/{name}_2.fastq",
                ]
            )
        elif layout == "SINGLE":
            cmds += cmd(
                [
                    "fastq-dump",
                    f"--outdir {out_dir}",
                    f"{srr_id}",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{out_dir}/{srr_id}.fastq",
                    f"{out_dir}/{name}.fastq",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{out_dir}/{name}.fastq",
                ]
            )

        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    if not check:
        df = pd.read_csv(sra_file, low_memory=False)
        tsc.logg.info(f"There are {df.shape[0]} samples.")

        if df["Library Name"].nunique() != df.shape[0]:
            tsc.logg.error(
                "Number of unique `Library Name` and total samples are not equal!"
            )

        df["cmd"] = df.apply(
            lambda x: cmds(x["Run"], x["Library Name"], x["LibraryLayout"]), axis=1
        )
        jobname = os.path.basename(__file__)[:-3]
        cmdmain = write_cmds_get_main(
            df,
            jobname,
            "12:00:00",
            "10",
            "python/3.7,sratoolkit/2.10.8",
            1,
            "farid.rsh@gmail.com",
            f"{out_dir}/_tmp",
        )
        os.system(cmdmain)
    else:
        tsc.logg.info("checking.")

    return None
