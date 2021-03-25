from pathlib import Path

from setuptools import find_packages, setup

try:
    from trisicell import __author__, __email__, __maintainer__, __version__
except ImportError:
    __author__ = ", ".join(["Farid Rashidi"])
    __maintainer__ = ", ".join(["Farid Rashidi"])
    __email__ = ", ".join(["farid.rsh@gmail.com"])
    __version__ = "0.0.1"

if __name__ == "__main__":
    setup(
        name="trisicell",
        entry_points="""
            [console_scripts]
            trisicell=trisicell.commands.trisicell:cli
        """,
        # ext_modules=extensions,
        # include_dirs=[np.get_include()],
        # package_data=package_data,
        use_scm_version=True,
        setup_requires=["setuptools_scm"],
        python_requires=">=3.6",
        install_requires=[
            l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
        ],
        extras_require=dict(
            dev=["pre-commit>=2.9.0"],
            docs=[
                l.strip()
                for l in (Path("docs") / "requirements.txt")
                .read_text("utf-8")
                .splitlines()
                if not l.startswith("-r")
            ],
        ),
        platforms=["Linux", "MacOSX"],
        packages=find_packages(),
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
