Installation
------------

trisicell requires Python 3.6 or later. We recommend to use Miniconda_.

PyPI
^^^^

Install trisicell from PyPI_ using::

    pip install -U trisicell

``-U`` is short for ``--upgrade``.
If you get a ``Permission denied`` error, use ``pip install -U trisicell --user`` instead.


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    pip install git+https://github.com/faridrashidi/trisicell

or::

    git clone https://github.com/faridrashidi/trisicell
    pip install -e trisicell

``-e`` is short for ``--editable`` and links the package to the original cloned
location such that pulled changes are also reflected in the environment.

To contribute to trisicell, ``cd`` into the cloned directory and
install the latest packages required for development together with the pre-commit hooks::

    pip install -r requirements-dev.txt
    pre-commit install


Dependencies
^^^^^^^^^^^^

- `numpy <https://numpy.org>`_, `scipy <https://scipy.org/>`_, `scikit-learn <https://scikit-learn.org/>`_.
- `pandas <https://pandas.pydata.org/>`_, `anndata <https://anndata.readthedocs.io/>`_.
- `matplotlib <https://matplotlib.org/>`_, `seaborn <https://seaborn.pydata.org>`_.


If you run into issues, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/trisicell
.. _Github: https://github.com/faridrashidi/trisicell
.. _`Github issue`: https://github.com/faridrashidi/trisicell/issues/new/choose
