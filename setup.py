import os
import sys
from pathlib import Path
from sys import platform

from setuptools import find_packages, setup
from setuptools.extension import Extension

try:
    from trisicell import __author__, __email__, __maintainer__, __version__
except ImportError:
    __author__ = ", ".join(["Farid Rashidi"])
    __maintainer__ = ", ".join(["Farid Rashidi"])
    __email__ = ", ".join(["farid.rsh@gmail.com"])
    __version__ = "0.0.13"

if platform == "linux" or platform == "linux2":
    os.environ["CC"] = "g++"
elif platform == "darwin":
    os.environ["CC"] = "clang++"
extensions = [
    Extension(
        "trisicell.external._mltd",
        sources=["trisicell/external/_mltd.pyx", "trisicell/external/mltd/mltd.cpp"],
        include_dirs=["trisicell/external/mltd"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "trisicell.external._scprob",
        sources=[
            "trisicell/external/_scprob.pyx",
            "trisicell/external/scprob/main.cpp",
        ],
        include_dirs=["trisicell/external/scprob"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "trisicell.external._scistree",
        sources=[
            "trisicell/external/_scistree.pyx",
            "trisicell/external/scistree/Utils.cpp",
            "trisicell/external/scistree/Utils2.cpp",
            "trisicell/external/scistree/Utils3.cpp",
            "trisicell/external/scistree/Utils4.cpp",
            "trisicell/external/scistree/UtilsNumerical.cpp",
            "trisicell/external/scistree/RerootTreeUtils.cpp",
            "trisicell/external/scistree/TreeBuilder.cpp",
            "trisicell/external/scistree/UnWeightedGraph.cpp",
            "trisicell/external/scistree/MarginalTree.cpp",
            "trisicell/external/scistree/RBT.cpp",
            "trisicell/external/scistree/PhylogenyTreeBasic.cpp",
            "trisicell/external/scistree/PhylogenyTree.cpp",
            "trisicell/external/scistree/BioSequenceMatrix.cpp",
            "trisicell/external/scistree/BinaryMatrix.cpp",
            "trisicell/external/scistree/GenotypeMatrix.cpp",
            "trisicell/external/scistree/ScistGenotype.cpp",
            "trisicell/external/scistree/ScistPerfPhyUtils.cpp",
            "trisicell/external/scistree/ScistPerfPhyImp.cpp",
            "trisicell/external/scistree/ScistDoublet.cpp",
            "trisicell/external/scistree/ScistErrRateInf.cpp",
            "trisicell/external/scistree/main.cpp",
        ],
        include_dirs=["trisicell/external/scistree"],
        extra_compile_args=["-O3", "-std=c++11", "-c"],
        language="c++",
    ),
    Extension(
        "trisicell.external._scite",
        sources=[
            "trisicell/external/_scite.pyx",
            "trisicell/external/scite/matrices.cpp",
            "trisicell/external/scite/mcmcBinTreeMove.cpp",
            "trisicell/external/scite/mcmc.cpp",
            "trisicell/external/scite/mcmcTreeMove.cpp",
            "trisicell/external/scite/output.cpp",
            "trisicell/external/scite/rand.cpp",
            "trisicell/external/scite/scoreBinTree.cpp",
            "trisicell/external/scite/scoreTree.cpp",
            "trisicell/external/scite/treelist.cpp",
            "trisicell/external/scite/trees.cpp",
            "trisicell/external/scite/findBestTrees.cpp",
        ],
        include_dirs=["trisicell/external/scite"],
        extra_compile_args=["-O3", "-std=c++11", "-c"],
        language="c++",
    ),
]


def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                sfile = path + ".cpp"
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0)))
if CYTHONIZE:
    try:
        from Cython.Build import cythonize
    except ImportError:
        sys.stderr.write(
            "Cannot find Cython. Have you installed all the requirements?\n"
            "Try pip install -r requirements.txt\n"
        )
        sys.exit(1)
    compiler_directives = {"language_level": 2, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)


setup(
    name="trisicell",
    entry_points="""
        [console_scripts]
        trisicell=trisicell.commands.trisicell:cli
    """,
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    python_requires=">=3.6",
    install_requires=[
        r.strip() for r in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    ext_modules=extensions,
    dependency_links=["https://pypi.gurobi.com"],
    extras_require={
        "dev": [
            "black==20.8b1",
            "pre-commit==2.9.3",
            "isort>=5.7.0",
            "pytest-cov",
        ],
        "docs": [
            r.strip()
            for r in (Path("docs") / "requirements.txt").read_text("utf-8").splitlines()
            if not r.startswith("-r")
        ],
    },
    platforms=["Linux", "MacOSX"],
    packages=find_packages(),
    include_package_data=True,
    author=__author__,
    author_email=__email__,
    email=__email__,
    maintainer=__maintainer__,
    maintainer_email=__email__,
    version=__version__,
    description=Path("README.rst").read_text("utf-8").split("\n")[3],
    long_description=Path("README.rst").read_text("utf-8"),
    long_description_content_type="text/x-rst; charset=UTF-8",
    license="BSD",
    url="https://github.com/faridrashidi/trisicell",
    project_urls={
        "Documentation": "https://trisicell.readthedocs.io/en/latest",
        "Source Code": "https://github.com/faridrashidi/trisicell",
    },
    download_url="https://github.com/faridrashidi/trisicell",
    keywords=[
        "tumor phylogeny",
        "single cell",
        "scalable",
        "rna-seq",
        "dna-seq",
    ],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Framework :: Jupyter",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
