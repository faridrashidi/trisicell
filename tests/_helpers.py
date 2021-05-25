import pytest

import trisicell as tsc

skip_gurobi = pytest.mark.skipif(
    tsc.ul.import_gurobi()[1], reason="Unable to import `Gurobi`!"
)
skip_mpi4py = pytest.mark.skipif(
    tsc.ul.import_mpi4py()[1], reason="Unable to import `MPI`!"
)
skip_rpy2 = pytest.mark.skipif(
    tsc.ul.import_rpy2()[1], reason="Unable to import `rpy2`!"
)
skip_graphviz = pytest.mark.skipif(
    tsc.ul.import_graphviz()[1], reason="Unable to import `Graphviz`!"
)
skip_graph_tool = pytest.mark.skipif(
    tsc.ul.import_graph_tool()[1], reason="Unable to import `graph_tool`!"
)
