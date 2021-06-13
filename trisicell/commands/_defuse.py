import glob

import click
import pandas as pd

from trisicell.ul._servers import *


@click.command(short_help="Run deFuse.")
@click.argument(
    "fq_files",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def defuse(fq_files):
    def cmds(subline):
        cmds = ""
        cmds += cmd([f"module load defuse/0.8.1"])
        cmds += cmd([f"mkdir -p /data/rashidimehrabf2/Fusion/map/{subline}"])
        cmds += cmd(
            [
                f"defuse.pl",
                f"-c /home/rashidimehrabf2/config_mm10_ens84.txt",
                f"-o /data/rashidimehrabf2/Fusion/map/{subline}",
                f"-d /fdb/defuse/mm10_ens84_newgmap",
                f"-1 {fq_files}/{subline}_R1_001.fastq",
                f"-2 {fq_files}/{subline}_R2_001.fastq",
                f"-s direct",
                f"-p $SLURM_CPUS_PER_TASK",
            ]
        )
        cmds += cmd([f"echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    temp = []
    for x in glob.glob(f"{fq_files}/*_R1_001.fastq"):
        temp.append(x.split("/")[-1].replace("_R1_001.fastq", ""))
    df["subline"] = temp
    df["cmd"] = df.apply(lambda x: cmds(x["subline"]), axis=1)

    cmdmain = write_cmds_get_main(
        df,
        tmpdir="/data/rashidimehrabf2/Fusion/tmp",
        jobname="rundefuse",
        time="2-00:00:00",
        mem="50",
        module="defuse/0.8.1",
        nthread=16,
    )
    os.system(cmdmain)

    return None
