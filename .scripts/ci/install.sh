#!/usr/bin/env bash

set -ev

echo "Installing APT dependencies"
if [[ "$OS" == "macos-latest" ]]; then
  brew install gcc
  brew install R
  brew install graphviz
  pip install pygraphviz
  brew install graph-tool

  python -m pip install --upgrade pip
  pip install -e .
  pip install pytest pytest-cov
elif [[ "$OS" == "ubuntu-latest" ]]; then
  python -m pip install --upgrade pip
  pip install -e .
  pip install pytest pytest-cov
  # https://github.com/yarnpkg/yarn/issues/7866
  curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | sudo apt-key add -
  # R-related
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
  sudo apt-get update -y
  sudo apt-get install libopenblas-base r-base r-base-dev -y
  sudo apt-get install libcurl4-openssl-dev libssl-dev -y
  pip install rpy2>=3.3.0
  sudo Rscript -e 'install.packages("devtools")'
  sudo Rscript -e 'install.packages("BiocManager")'
  sudo Rscript -e 'install.packages("tidyverse")'
  sudo Rscript -e 'install.packages("cowplot")'
  sudo Rscript -e 'devtools::install_github("YuLab-SMU/ggtree")'
  sudo Rscript -e 'devtools::install_github("YuLab-SMU/aplot")'
  sudo Rscript -e 'devtools::install_github("xiangpin/ggtreeExtra")'
  sudo Rscript -e 'devtools::install_github("zhouzilu/DENDRO")'
  sudo Rscript -e 'BiocManager::install("graph")'
  sudo Rscript -e 'devtools::install_bitbucket("edith_ross/oncoNEM")'

  sudo apt-get install graphviz graphviz-dev -y
  pip install pygraphviz

  sudo add-apt-repository "deb [ arch=amd64 ] https://downloads.skewed.de/apt bullseye main"
  sudo apt-key adv --keyserver keys.openpgp.org --recv-key 612DEFB798507F25
  sudo apt-get install python3-graph-tool -y
else
  exit 42
fi