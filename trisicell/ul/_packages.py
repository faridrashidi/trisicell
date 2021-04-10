import trisicell as tsc


def import_gurobi():
    try:
        import gurobipy as gp

        # m = gp.Model("MIP")
        # m.Params.OutputFlag = 0
        # m.Params.LogFile = ""
        # m.Params.Threads = 1
        # x = m.addVar(vtype=gp.GRB.BINARY, name="x")
        # y = m.addVar(vtype=gp.GRB.BINARY, name="y")
        # z = m.addVar(vtype=gp.GRB.BINARY, name="z")
        # m.setObjective(x + y + 2 * z, gp.GRB.MAXIMIZE)
        # m.addConstr(x + 2 * y + 3 * z <= 4, "c0")
        # m.addConstr(x + y >= 1, "c1")
        # m.optimize()
        return gp, False
    except:
        tsc.logg.error(
            "Unable to import `gurobipy`!",
            "Make sure `Gurobi` is already installed in your system.",
            "Then install `gurobipy`.",
        )
        return None, True


def import_mpi4py():
    try:
        import mpi4py

        return mpi4py, False
    except:
        tsc.logg.error(
            "Unable to import `mpi4py`!",
            "Make sure `mpi/mpich-x86_64` is already installed in your system.",
            "Then install `mpi4py` using `pip`.",
        )
        return None, True


def import_rpy2(name="base", how=""):
    try:
        from rpy2.robjects import r
        from rpy2.robjects.packages import PackageNotInstalledError, importr

    except ImportError:
        tsc.logg.error(
            "Unable to import `rpy2`, install it first as `pip install rpy2` version `>=3.3.0`."
        )
        return None, True

    try:
        _r_lib = importr(name)
        return _r_lib, False
    except PackageNotInstalledError:
        tsc.logg.error(f"Install R library `{name!r}` first." f"{how}")
        return None, True


def import_graphviz():
    try:
        import pygraphviz

        return pygraphviz, False
    except:
        tsc.logg.error(
            "Unable to import `pygraphviz`!",
            "Make sure `Graphviz` is already installed in your system.",
            "Then install `pygraphviz` using `pip`.",
        )
        return None, True
