CLI
===

A command line interface (CLI) is available in Trisicell package.
After you have trisicell correctly installed on your machine
(see :ref:`installation tutorial <installationguide>`), the ``trisicell``
command will become available in the terminal. ``trisicell`` is a
command line tool with subcomands. You can get quick info on all the
available commands typing ``trisicell --help``. You will get the
following output:

.. code-block:: bash

    Usage: trisicell [OPTIONS] COMMAND [ARGS]...

      Trisicell.

      Scalable intratumor heterogeneity inference and validation from single-cell data.

    Options:
      --version  Show the version and exit.
      --help     Show this message and exit.

    Commands:
      mcalling   Mutation calling.
      booster    Boost available tree reconstruction tool (Trisicell-Boost).
      partf      Get samples or calculate for PartF.
      consensus  Build consensus tree between two phylogenetic trees (Trisicell-Cons).


``mcalling`` - Run Mutation Calling
-----------------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: mcalling
    :nested: full


``booster`` - Run Booster
-------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: booster
    :nested: full


``consensus`` - Run Consensus
-----------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: consensus
    :nested: full
