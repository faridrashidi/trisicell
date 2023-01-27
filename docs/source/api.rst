.. module:: trisicell
.. automodule:: trisicell
    :noindex:

API
===

Import Trisicell as::

   import trisicell as tsc

After mutation calling and building the input data via our suggested
:ref:`mutation calling pipeline <caller>`.


Read/Write (io)
---------------
This module offers a bunch of functions for reading and writing of the data.

.. module:: trisicell.io
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    io.read
    io.write


Preprocessing (pp)
------------------
This module offers a bunch of functions for filtering and preprocessing of the
data.

.. module:: trisicell.pp
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    pp.remove_mut_by_list
    pp.remove_cell_by_list
    pp.filter_mut_reference_must_present_in_at_least
    pp.filter_mut_mutant_must_present_in_at_least
    pp.consensus_combine


Tools (tl)
----------
This module offers a high-level API to compute the conflict-free solution
and calculating the probability of mutations seeding particular cells.

.. module:: trisicell.tl
.. currentmodule:: trisicell

**Solving the noisy input genotype matrix (Trisicell-Boost)**

.. autosummary::
    :toctree: .

    tl.booster

**Partition function calculation (Trisicell-PartF)**

.. autosummary::
    :toctree: .

    tl.partition_function

**Consensus tree building (Trisicell-Cons)**

.. autosummary::
    :toctree: .

    tl.consensus


Plotting (pl)
-------------
This module offers plotting the tree in clonal or dendrogram format.

.. module:: trisicell.pl
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    pl.clonal_tree
    pl.dendro_tree


Utils (ul)
----------
This module offers a bunch of utility functions.

.. module:: trisicell.ul
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    ul.to_tree
    ul.to_cfmatrix
    ul.to_mtree


Datasets (datasets)
-------------------
This module offers a bunch of functions for simulating data.

.. module:: trisicell.datasets
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    datasets.example
    datasets.sublines_bwes
    datasets.sublines_bwts
    datasets.sublines_scrnaseq
