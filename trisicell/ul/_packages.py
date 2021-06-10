import trisicell as tsc


def import_gurobi():
    try:
        import gurobipy as gp

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
            "Unable to import `rpy2`, install it first as `pip install rpy2` version"
            " `>=3.3.0`."
        )
        return None, True

    try:
        _r_lib = importr(name)
        return _r_lib, False
    except PackageNotInstalledError:
        tsc.logg.error(f"Install R library `{name!r}` first.{how}")
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


def import_graph_tool():
    try:
        import graph_tool

        return graph_tool, False
    except:
        tsc.logg.error(
            "Unable to import `graph_tool`!",
            "Make sure `graph_tool` is installed in your system.",
        )
        return None, True
