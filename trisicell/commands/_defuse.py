import glob
import os

import click
import pandas as pd

import trisicell as tsc
from trisicell.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run deFuse.")
@click.argument(
    "in_files",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "out_files",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
def defuse(in_files, out_files):
    """Run deFuse.

    trisicell defuse PATH_TO_IN_DIR PATH_TO_OUT_DIR
    """

    def cmds(cell):
        cmds = ""
        cmds += cmd(["module load defuse/0.8.1"])
        cmds += cmd([f"mkdir -p {out_files}/{cell}"])
        cmds += cmd(
            [
                "defuse.pl",
                "-c /home/rashidimehrabf2/config_mm10_ens84.txt",
                f"-o {out_files}/{cell}",
                "-d /fdb/defuse/mm10_ens84_newgmap",
                f"-1 {in_files}/{cell}_R1_001.fastq",
                f"-2 {in_files}/{cell}_R2_001.fastq",
                "-s direct",
                "-p $SLURM_CPUS_PER_TASK",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    temp = []
    for x in glob.glob(f"{in_files}/*_R1_001.fastq"):
        temp.append(x.split("/")[-1].replace("_R1_001.fastq", ""))
    df["cell"] = temp
    df["cmd"] = df.apply(lambda x: cmds(x["cell"]), axis=1)
    tsc.logg.info(f"Number of cells: {df.shape[0]}")

    jobname = os.path.basename(__file__)[:-3]
    cmdmain = write_cmds_get_main(
        df,
        jobname,
        "2-00:00:00",
        "50",
        "defuse/0.8.1",
        16,
        "farid.rsh@gmail.com",
        f"{out_files}/_tmp",
    )
    os.system(cmdmain)

    return None
