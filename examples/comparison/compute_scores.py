"""
Comparing scores for two phylogenetic trees
-------------------------------------------

This example shows how to compare/measure two inferred genotype data (trees).
"""

import trisicell as tsc

# %%
# First, we load two binary test single-cell genotype data.
grnd = tsc.io.read(
    tsc.ul.get_file("trisicell.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
)
sol = tsc.io.read(
    tsc.ul.get_file("trisicell.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
)

# %%
# Calculating the ancestor-descendent accuracy.
tsc.tl.ad(grnd, sol)

# %%
# Calculating the different-lineage accuracy.
tsc.tl.dl(grnd, sol)

# %%
# Calculating the multi-labeled tree dissimilarity measure (MLTD).
tsc.tl.mltd(grnd, sol)

# %%
# Calculating the tumor phylogeny tree edit distance measure (TPTED).
tsc.tl.tpted(grnd, sol)

# %%
# Calculating the distinctly inherited sets score (DISC).
tsc.tl.disc(grnd, sol)

# %%
# Calculating the commonly ancestor sets score (CASet).
tsc.tl.caset(grnd, sol)

# %%
# Calculating the Triplet-based similarity score (MP3).
tsc.tl.mp3(grnd, sol)
