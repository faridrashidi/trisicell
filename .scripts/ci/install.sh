#!/usr/bin/env bash

set -ev

python -m pip install --upgrade pip
pip install -e .
pip install pytest pytest-cov

if [[ "$OS" == "macos-latest" ]]; then
  brew install R
  brew install graphviz
  pip install pygraphviz
  brew install graph-tool
elif [[ "$OS" == "ubuntu-latest" ]]; then
  echo "Installing APT dependencies"
  # https://github.com/yarnpkg/yarn/issues/7866
  curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | sudo apt-key add -
  # R-related
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
  sudo apt-get update -y
  sudo apt-get install libopenblas-base r-base r-base-dev -y
  pip install rpy2>=3.3.0
  Rscript -e 'install.packages("devtools")'
  Rscript -e 'install.packages("BiocManager")'
  Rscript -e 'install.packages("tidyverse")'
  Rscript -e 'install.packages("cowplot")'
  Rscript -e 'devtools::install_github("YuLab-SMU/ggtree")'
  Rscript -e 'devtools::install_github("YuLab-SMU/aplot")'
  Rscript -e 'devtools::install_github("xiangpin/ggtreeExtra")'
  Rscript -e 'devtools::install_github("zhouzilu/DENDRO")'
  Rscript -e 'BiocManager::install("graph")'
  Rscript -e 'devtools::install_bitbucket("edith_ross/oncoNEM")'

  sudo apt-get install graphviz graphviz-dev -y
  pip install pygraphviz
  # Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi