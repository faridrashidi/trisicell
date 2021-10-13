.. module:: trisicell
.. automodule:: trisicell
    :noindex:

API
===

Import Trisicell as::

   import trisicell as tsc

After mutation calling and building the input data via our suggested
:ref:`mutation calling pipeline <caller>`.


Datasets (datasets)
-------------------
This module offers a bunch of functions for simulating data.

.. module:: trisicell.datasets
.. currentmodule:: trisicell
.. autosummary::
    :toctree: .

    datasets.example
    datasets.simulate
    datasets.add_noise
    datasets.melanoma20
    datasets.colorectal1
    datasets.colorectal2
    datasets.colorectal3
    datasets.acute_lymphocytic_leukemia1
    datasets.acute_lymphocytic_leukemia2
    datasets.acute_lymphocytic_leukemia3
    datasets.acute_lymphocytic_leukemia4
    datasets.acute_lymphocytic_leukemia5
    datasets.acute_lymphocytic_leukemia6
    datasets.high_grade_serous_ovarian_cancer1
    datasets.high_grade_serous_ovarian_cancer2
    datasets.high_grade_serous_ovarian_cancer3
    datasets.high_grade_serous_ovarian_cancer_3celllines
    datasets.myeloproliferative_neoplasms18
    datasets.myeloproliferative_neoplasms78
    datasets.myeloproliferative_neoplasms712
    datasets.renal_cell_carcinoma
    datasets.muscle_invasive_bladder
    datasets.erbc
    datasets.tnbc


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
    pp.bifiltering
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
    tl.scite
    tl.phiscsb
    tl.scistree
    tl.onconem
    tl.huntress


**Partition function calculation (Trisicell-PartF)**

.. autosummary::
    :toctree: .

    tl.partition_function

**Consensus tree building (Trisicell-Cons)**

.. autosummary::
    :toctree: .

    tl.consensus

**For comparing two phylogenetic trees**

.. autosummary::
    :toctree: .

    tl.ad
    tl.dl
    tl.mltd
    tl.tpted
    tl.caset
    tl.disc
    tl.mp3
    tl.rf


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
    ul.hclustering
    ul.is_conflict_free_gusfield
