#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021, Farid Rashidi Mehrabadi All rights reserved.

# =========================================================================================
# Author     : Farid Rashidi Mehrabadi (farid.rashidimehrabadi@nih.gov)
# Last Update: Aug 14, 2020
# Description: cleaning
# =========================================================================================

import glob
import subprocess

from trisicell.ul._servers import *


def _is_ok(name):
    file = open(name, "r")
    body = file.read()
    file.close()
    a = body.count("&& echo Done! )")
    b = body.count("Done!\n")
    if a == 0 and b == 1:
        return True
    else:
        return a == b


def after01(config):
    if config["isrna"] == True:
        steps = [
            "s01indexing",
            "s02mapping",
            "s03indexing",
            "s04mapping",
            "s05calling",
            "s06jointcalling",
            "s07merging",
            "s08annotating",
            "s09expressing",
            "s10velocitying",
        ]
    else:
        steps = [
            "s02mapping",
            "s04mapping",
            "s05calling",
            "s06jointcalling",
            "s07merging",
            "s08annotating",
        ]

    conds = {}
    for cond in steps:
        x = 0
        for file in glob.glob(f"{config['tmpdir']}/log/{cond}/*.o"):
            if not _is_ok(file):
                x += 1
        conds[cond] = x
    print(conds)
