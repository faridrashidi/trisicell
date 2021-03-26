CLI
===

A command line interface (CLI) is available in Trisicell package. After you have trisicell correctly
installed on your machine (see :ref:`installation tutorial <installation>`), the ``trisicell`` command will
become available in the terminal.
``trisicell`` is a command line tool with subcomands. You can get quick info on all the available commands
typing ``trisicell --help``. You will get the following output:

.. code-block:: bash
    
    Usage: trisicell [OPTIONS] COMMAND [ARGS]...

        Scalable tumor phylogeny reconstruction and validation

        Options:
            --version  Show the version and exit.
            --help     Show this message and exit.

        Commands:
            score     Caculate scores.
            scistree  Run ScisTree.
            scite     Run SCITE.
            booster   Run Booster.
            phiscsb   Run PhISCS (CSP version).
            phiscsi   Run PhISCS (ILP version).
            bnb       Run PhISCS-BnB.
            huntress  Run HUNTRESS.

``scistree`` - Run ScisTree
---------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: scistree
    :nested: full


``scite`` - Run SCITE
---------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: scite
    :nested: full


``score`` - Calculating Scores
------------------------------

.. click:: trisicell.commands.trisicell:cli
    :prog: trisicell
    :commands: score
    :nested: full
