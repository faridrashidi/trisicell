import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import cut_tree, dendrogram, fcluster, linkage


def calc_hclust(dist, names, method="ward"):
    clust = linkage(dist[np.triu_indices(dist.shape[0], 1)], method=method)
    clusters = {}

    for i in range(2, dist.shape[0]):
        fc = fcluster(clust, i, criterion="maxclust")
        clusters[i] = pd.Series(fc, index=names)

    return clusters

    # sns.set(style='ticks', font='Helvetica', font_scale=0.4)
    # dend = dendrogram(clust, color_threshold=threshold, labels=names)

    # dend = dendrogram(clust)
    # idx = [int(x) for x in dend["ivl"]]
    # corrcoef = 1 - dist
    # sns.clustermap(
    #     corrcoef[idx, :][:, idx],
    #     vmin=vmin,
    #     vmax=vmax,
    #     cmap="RdYlBu_r",
    #     row_cluster=False,
    #     col_cluster=False,
    #     cbar_pos=(1, 0.03, 0.05, 0.94),
    #     figsize=(5, 5),
    #     xticklabels=names[idx],
    #     yticklabels=False,
    #     dendrogram_ratio=0.000001,
    # )
