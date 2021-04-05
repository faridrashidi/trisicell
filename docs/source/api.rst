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
~~~~~~~~~~~~~~~
This module offers a bunch of functions for reading and writing of the data.

.. module:: trisicell.io
.. currentmodule:: trisicell
.. autosummary::
    :toctree: gen_modules

    io.read
    io.read_gatk
    io.read_rsem
    io.write


Datasets (datasets)
~~~~~~~~~~~~~~~~~~~
This module offers a bunch of functions for simulating data.

.. module:: trisicell.datasets
.. currentmodule:: trisicell
.. autosummary::
    :toctree: gen_modules

    datasets.example
    datasets.simulate
    datasets.add_noise
    datasets.melanoma20
    datasets.colorectal2
    datasets.acute_lymphocytic_leukemia2
    datasets.high_grade_serous_ovarian_cancer


Preprocessing (pp)
~~~~~~~~~~~~~~~~~~
This module offers a bunch of functions for filtering and preprocessing of the data.

.. module:: trisicell.pp
.. currentmodule:: trisicell
.. autosummary::
    :toctree: gen_modules

    pp.remove_mut_by_list
    pp.remove_cell_by_list
    pp.filter_mut_reference_must_present_in_at_least
    pp.filter_mut_mutant_must_present_in_at_least


Tools (tl)
~~~~~~~~~~
This module offers a high-level API to compute the conflict-free solution
and calculating the probability of mutations seeding particular cells.

.. module:: trisicell.tl
.. currentmodule:: trisicell

Solving the noisy input genotype matrix (Trisicell-Boost)
#########################################################

.. autosummary::
    :toctree: gen_modules

    tl.solver.booster
    tl.solver.scite
    tl.solver.phiscsb
    tl.solver.scistree
    

Partition function calculation (Trisicell-PartF)
################################################

.. autosummary::
    :toctree: gen_modules

    tl.partition_function

Consensus tree building (Trisicell-Cons)
########################################

.. autosummary::
    :toctree: gen_modules

    tl.consensus_combine
    tl.consensus_tree

For comparing two phylogenetic trees
########################################

.. autosummary::
    :toctree: gen_modules

    tl.score.ad
    tl.score.dl
    tl.score.mltd
    tl.score.tpted


Plotting (pl)
~~~~~~~~~~~~~
This module offers plotting the tree in clonal or dendrogram format.

.. module:: trisicell.pl
.. currentmodule:: trisicell
.. autosummary::
    :toctree: gen_modules

    pl.clonal_tree
    pl.dendro_tree


Utils (ul)
~~~~~~~~~~
This module offers a bunch of utility functions.

.. module:: trisicell.ul
.. currentmodule:: trisicell
.. autosummary::
    :toctree: gen_modules

    ul.is_conflict_free_gusfield
    ul.to_tree
