.. automodule:: trisicell

API
===

Import Trisicell as::

   import trisicell as tsc

After reading the data...


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

**Solving the noisy input genotype matrix**

.. autosummary::
    :toctree: gen_modules

    tl.solver.dnc
    tl.solver.scistree
    tl.solver.scite
    tl.solver.phiscs
    tl.solver.cardelino

**Partition function calculation (Trisicell-PF)**

.. autosummary::
    :toctree: gen_modules

    tl.partf


Plotting (pl)
~~~~~~~~~~~~~
This module offers plotting the tree in dendrogram or clonal format.

.. module:: trisicell.pl
.. currentmodule:: trisicell

.. autosummary::
    :toctree: gen_modules

    pl.get_tree


Utils (ul)
~~~~~~~~~~
This module offers a bunch of utility functions.

.. module:: trisicell.ul
.. currentmodule:: trisicell

.. autosummary::
    :toctree: gen_modules

    ul.is_conflict_free_gusfield
    ul.count_flips
    ul.infer_rates
