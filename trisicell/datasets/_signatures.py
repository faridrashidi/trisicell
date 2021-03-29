import pandas as pd
import yaml

import trisicell as tsc


def get_markers():
    markers = pd.read_csv(tsc.ul.get_file("trisicell.datasets/data/markers.csv"))
    return markers


def get_signatures(kind):
    signature_file = ""
    if kind == "mm10":
        signature_file = tsc.ul.get_file("trisicell.datasets/data/signatures_mm10.yml")
    elif kind == "hg19":
        signature_file = tsc.ul.get_file("trisicell.datasets/data/signatures_hg19.yml")
    with open(signature_file) as fin:
        signatures = yaml.load(fin, Loader=yaml.FullLoader)
    return signatures
