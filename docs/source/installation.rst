:tocdepth: 1

.. _installationguide:

Installation
------------

Trisicell requires Python 3.6 or later.

Using PyPI
^^^^^^^^^^

Install trisicell from PyPI_ using::

    pip install -U trisicell

``-U`` is short for ``--upgrade``.
If you get a ``Permission denied`` error, use
``pip install -U trisicell --user`` instead.


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    git clone https://github.com/faridrashidi/trisicell
    cd trisicell
    git fetch
    git checkout develop
    pip install -e '.[dev]'

``-e`` stands for ``--editable`` and makes sure that your environment is
updated when you pull new changes from GitHub. The ``'[dev]'`` options
installs requirements needed for development, because Trisicell is bundled
with an additional library.

Your contributions to improve Trisicell is highly appreciated! Please
check out our `contributing guide`_.

If you run into issues, do not hesitate to approach us or raise a
`GitHub issue`_.

.. _PyPI: https://pypi.org/project/trisicell
.. _Github: https://github.com/faridrashidi/trisicell
.. _`Github issue`: https://github.com/faridrashidi/trisicell/issues/new/choose
.. _contributing guide: https://github.com/faridrashidi/trisicell/blob/master/CONTRIBUTING.rst
