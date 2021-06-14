# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from datetime import datetime
from pathlib import Path

from pybtex.plugin import register_plugin
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.labels import BaseLabelStyle
from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey

import trisicell

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, os.path.abspath("_ext"))

# -- Retrieve notebooks ------------------------------------------------
# for nb in glob(os.path.join(".", "*.ipynb")):
#     os.remove(nb)
# notebooks_url = "../../trisicell_notebooks/"
# notebooks = [
#     "Mel.Subline.ipynb",
#     "Mel.Treated.ipynb",
# ]
# for nb in notebooks:
#     copyfile(os.path.join(notebooks_url, nb), os.path.join(".", nb))


# -- Project information -----------------------------------------------------
nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = "2.0"  # Nicer param docs
suppress_warnings = ["ref.citation"]
project = "Trisicell"
author = "National Cancer Institute"  # trisicell.__author__
title = (
    "Scalable intratumor heterogeneity inference and validation from single-cell data"
)
copyright = f"{datetime.now():%Y}, {author}"
release = "master"
version = f"master ({trisicell.__version__})"


# -- General configuration ---------------------------------------------------
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    "sphinx_autodoc_typehints",
    "nbsphinx",
    "sphinx_paramlinks",
    "sphinx_copybutton",
    "sphinx_click.ext",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx_gallery.gen_gallery",
    # "sphinx_last_updated_by_git",
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False

# bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "mystyle"

# spelling
spelling_lang = "en_US"
spelling_warning = True
spelling_add_pypi_package_names = True
spelling_show_suggestions = True
spelling_exclude_patterns = ["references.rst"]
# see: https://pyenchant.github.io/pyenchant/api/enchant.tokenize.html
spelling_filters = [
    "enchant.tokenize.URLFilter",
    "enchant.tokenize.EmailFilter",
    "docs.source.utils.ModnameFilter",
    "docs.source.utils.SignatureFilter",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = [".rst", ".ipynb"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
master_doc = "index"
default_role = "literal"
pygments_style = "sphinx"
todo_include_todos = False

intersphinx_mapping = dict(
    python=("https://docs.python.org/3", None),
    numpy=("https://docs.scipy.org/doc/numpy/", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    networkx=("https://networkx.github.io/documentation/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    matplotlib=("https://matplotlib.org/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    anndata=("https://anndata.readthedocs.io/en/latest/", None),
)


# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
github_repo = "trisicell"
github_nb_repo = "trisicell_notebooks"
html_logo = "_static/images/logo.svg"
html_theme_options = {
    "logo_only": True,
    "display_version": True,
    "navigation_depth": 2,
}
html_context = {
    "display_github": True,
    "github_user": "faridrashidi",
    "github_repo": "trisicell",
    "github_version": "master",
    "conf_py_path": "/docs/",
}
html_show_sphinx = False
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]


# -- Options for other output ------------------------------------------
htmlhelp_basename = f"{project}doc"
title_doc = f"{project} documentation"

latex_documents = [(master_doc, f"{project}.tex", title_doc, author, "manual")]
man_pages = [(master_doc, project, title_doc, [author], 1)]
texinfo_documents = [
    (master_doc, project, title_doc, author, project, title, "Miscellaneous")
]


# a simple label style which uses the bibtex keys for labels
class MyLabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key


class MyStyle(UnsrtStyle):

    default_label_style = MyLabelStyle


register_plugin("pybtex.style.formatting", "mystyle", MyStyle)


# -- sphinx gallery ------------------------------------------
def reset_matplotlib(gallery_conf, fname):
    import matplotlib as mpl

    mpl.use("agg")

    import matplotlib.pyplot as plt

    plt.rcdefaults()
    mpl.rcParams["savefig.bbox"] = "tight"
    mpl.rcParams["savefig.transparent"] = True


example_dir = HERE.parent.parent / "examples"
rel_example_dir = Path("..") / ".." / "examples"


sphinx_gallery_conf = {
    "image_scrapers": "matplotlib",
    "reset_modules": (
        "seaborn",
        reset_matplotlib,
    ),
    "filename_pattern": f"{os.path.sep}(plot_|compute_)",
    "examples_dirs": example_dir,
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
    "abort_on_example_error": True,
    "show_memory": True,
    "within_subsection_order": FileNameSortKey,
    "subsection_order": ExplicitOrder(
        [
            rel_example_dir / "reconstruction",  # really must be relative
            rel_example_dir / "comparison",
        ]
    ),
    "reference_url": {
        "sphinx_gallery": None,
    },
    "line_numbers": False,
    "compress_images": ("images", "thumbnails"),
    "inspect_global_variables": False,
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": "trisicell",
    "download_all_examples": False,
    "pypandoc": True,  # convert rST to md when downloading notebooks
}
