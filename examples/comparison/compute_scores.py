"""
Comparing scores for two phylogenetic trees
-------------------------------------------

This example shows how to calculate ancestor-descendant, different-lineage accuracies.
"""

import trisicell as tsc

# %%
# First, we load a binary test single-cell genotype data.
grnd = tsc.io.read(
    tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
)
sol = tsc.io.read(
    tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
)
ad = tsc.tl.ad(grnd, sol)
ad
