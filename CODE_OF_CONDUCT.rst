Contributing guidelines
~~~~~~~~~~~~~~~~~~~~~~~

Welcome to `Trisicell project <https://github.com/faridrashidi/trisicell>`_!
Before sending your pull requests, make sure that you **read the whole
guidelines**. If you have any doubt on the contributing guide, please
feel free to `state it clearly in an
issue <https://github.com/faridrashidi/trisicell/issues/new/choose>`_.

Table of Contents
=================
- `Contributing to Trisicell`_
- `Codebase structure`_
- `Code style guide`_
- `Testing`_
- `Writing documentation`_
- `Building on a new machine`_
- `Building cpp files`_
- `Submitting a PR`_
- `Creating a new release`_


Contributing to Trisicell
-------------------------
Clone Trisicell from source as::

    git clone https://github.com/faridrashidi/trisicell
    cd trisicell

Install the development mode::

    pip install -e '.[dev]'

Then install pre-commit. This will ensure that the pushed code passes the
linting steps::

    pre-commit install


Codebase structure
------------------
The Trisicell project:

- `trisicell <trisicell>`_: the root of the package.

    - `trisicell/io <trisicell/io>`_: the read/write module, offers a bunch of
        functions for reading and writing of the data.
    - `trisicell/pl <trisicell/pl>`_: the plotting module, offers plotting the
        tree in clonal or dendrogram format.
    - `trisicell/pp <trisicell/pp>`_: the preprocessing module, offers a bunch
        of functions for filtering and preprocessing of the data.
    - `trisicell/tl <trisicell/tl>`_: the tools module, offers a high-level API
        to compute the conflict-free solution and calculating the probability of
        mutations seeding particular cells.
    - `trisicell/ul <trisicell/ul>`_: the utils module, offers a bunch of
        utility functions.
    - `trisicell/commands <trisicell/commands>`_: the CLI commands module,
        offers running trisicell in command-line interface (CLI) mode.
    - `trisicell/datasets <trisicell/datasets>`_: the datasets module, offers
        some of the published single-cell datasets and generating simulations.

Tests structure:

- `tests <tests>`_: the root of the test files.


Code style guide
----------------
We rely on ``black`` and ``isort`` to do the most of the formatting - both of
them are integrated as pre-commit hooks. You can use ``pre-commit`` to check
the changes::

    make lint

Keep in mind to use tags like ``TODO[colon]`` and ``FIXME[colon]`` to make
other developers aware of what to do next.


Testing
-------
We use ``pytest`` to automate our testing. To execute the tests, run::

    make test


Writing documentation
---------------------
We use ``numpy``-style docstrings for the documentation with the following
additions and modifications:

- when referring to some argument within the same docstring, enclose that
    reference in \`\`.
- prefer putting references in the ``references.bib`` instead under the
    ``References`` sections of the docstring.

In order to build the documentation, run::

    make doc


Building on a new machine
-------------------------
For building the ``.so`` files you need to execute::

    pip install -r requirements.txt
    python setup.py build
    python setup.py build_ext --inplace
    pip install -e .


Building cpp files
------------------
For building ``.cpp`` files from ``.pyx`` files you need to execute::

    pip install -r requirements.txt
    CYTHONIZE=1 python setup.py install
    pip install -e .


Submitting a PR
---------------
Before submitting a new pull request, please make sure you followed these
instructions:

- make sure that your code follows the above specified conventions (see
    `Code style guide`_ and `Writing documentation`_).
- if applicable, make sure you've added/modified at least 1 test to account
    for the changes you've made.
- make sure that all tests pass locally (see `Testing`_).
- if there is no issue which this PR solves, create a new
    `one <https://github.com/faridrashidi/trisicell/issues/new>`_ and briefly
    explaining what the problem is.


Creating a new release
----------------------
If you are a core developer and you want to create a new release, you need to
install ``bump2version`` first as::

    pip install bump2version

Depending on what part of the release you want to update, you can run::

    bump2version {major,minor,patch}

By default, this will create a new tag and automatically update the
``__version__`` wherever necessary, commit the changes and create a new tag.
If you have uncommitted files in the tree, you can use ``--allow-dirty``
flag to include them in the commit.

After the version has been bumped, make sure to push the commit **AND**
the newly create tag to the upstream. This can be done by e.g. setting
``push.followtags=true`` in your git config or use
``git push --atomic <branch> <tag>``.
