#!/usr/bin/env bash

set -ev

python -m pip install --upgrade pip
pip install -e .
pip install pytest pytest-cov

if [[ "$OS" == "macos-latest" ]]; then
  brew install graphviz
elif [[ "$OS" == "ubuntu-latest" ]]; then
  pip install rpy2>=3.3.0
  # Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi