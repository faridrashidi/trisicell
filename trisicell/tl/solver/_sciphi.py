import os
import time

import trisicell as tsc


def sciphi(df_input):
    # TODO: implement
    executable = tsc.ul.executable("sciphi", "SCIPhI")

    tsc.logg.info("running SCIPhI with")

    # tmpdir = tsc.ul.tmpdirsys(suffix=".sciphi")
    tmpdir = "test"
    tsc.ul.cleanup(tmpdir)
    tsc.ul.mkdir(tmpdir)

    matrix_I = df_input.values
    with open(f"{tmpdir}/sciphi.mpileup", "w") as fout:
        for j in range(matrix_I.shape[1]):
            line = f"seq1\t{(j+1)*100}\tA"
            r = q = ""
            for i in range(matrix_I.shape[0]):
                if matrix_I[i, j] == 0:
                    r = "."
                elif matrix_I[i, j] == 1:
                    r = "T"
                elif matrix_I[i, j] == 3:
                    r = "N"
                q = "<"
                line = f"{line}\t1\t{r}\t{q}"
            fout.write(line + "\n")
    with open(f"{tmpdir}/sciphi.cellnames", "w") as fout:
        for i in range(matrix_I.shape[0]):
            fout.write(f"{df_input.index[i]}\tCT\n")

    cmd = (
        f"{executable} "
        f"-o {tmpdir}/out "
        f"--in {tmpdir}/sciphi.cellnames "
        "--seed 42 "
        f"{tmpdir}/sciphi.mpileup "
        f"> {tmpdir}/sciphi.log"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time
    running_time

    return None
