"""
Construct lienage tree using SCITE
----------------------------------

This example shows how to construct a lineage tree using SCITE on a binary single-cell genotype matrix.
"""

import trisicell as tsc

# %%
# First, we load a binary test single-cell genotype data.
df_in = tsc.datasets.test()
df_in.head()

# %%
# Next, using :func:`trisicell.tl.scite` we remove the single-cell noises from the input.
df_out = tsc.tl.scite(df_in, alpha=0.0000001, beta=0.1, n_restarts=3, n_iters=1000)
df_out.head()

# %%
# Finally, using :func:`trisicell.ul.is_conflict_free_gusfield` we check whether the inferred
# genotype matrix is conflict-free or not.
is_cf = tsc.ul.is_conflict_free_gusfield(df_out)
is_cf
