import matplotlib.colors as mcolors
import natsort as ns
import numpy as np
import seaborn as sns


def heatmap(adata, color_attrs=None, layer="X", figsize=(12, 7)):
    # TODO: check for test
    if color_attrs is not None:
        row_colors = []
        for attr in color_attrs:
            row_colors.append(adata.obs[attr])
    else:
        row_colors = None

    if layer == "X":
        rvb = mcolors.LinearSegmentedColormap.from_list(
            "",
            ["#FFFFFF", "#000000"],  # "#A6CEE3"
        )
        df = adata.to_df().copy()
        df[df == 3] = 0
        vmin = 0
        vmax = 1
        col_cluster = True
    else:
        rvb = "RdBu_r"
        df = adata.obsm[layer].copy()
        vmin = -0.5
        vmax = 0.5

        # df[df > 0.2] = 1
        # df[(-0.2 <= df) & (0.2 >= df)] = 0
        # df[df < -0.2] = -1
        # vmin = -2
        # vmax = 2

        col_cluster = False

    sns.clustermap(
        df,
        vmin=vmin,
        vmax=vmax,
        metric="euclidean",
        cmap=rvb,
        row_cluster=False,
        col_cluster=col_cluster,
        row_colors=row_colors,
        cbar_pos=None,
        figsize=figsize,
        xticklabels=False,
        yticklabels=False,
        # linecolor="white",
        # linewidths=0.8,
        # annot=False,
        # colors_ratio=0.05,
        dendrogram_ratio=0.0001,
    )
    # plt.savefig(filepath, bbox_inches="tight", pad_inches=0)


def plot_flips(df_in, df_out, row_colors):
    # TODO: check for test
    snvs = np.intersect1d(df_in.columns, df_out.columns)
    # cells = np.intersect1d(df_in.index, df_out.index)
    tmp = []
    for x in snvs:
        tmp.append(".".join(x.split(".chr")[1].split(".")[:2]))
    tmp = ns.index_natsorted(tmp)
    snvs = snvs[tmp]
    I_mtr = df_in.loc[row_colors.keys(), snvs]
    O_mtr = df_out.loc[row_colors.keys(), snvs]

    ccolors = []
    chrs = set()
    for x in snvs:
        c = x.split(".chr")[1].split(".")[0]
        if c in [f"{i}" for i in range(1, 22, 2)] + ["Y"]:
            ccolors.append("#969696")
            chrs.add(c)
        if c in [f"{i}" for i in range(2, 22, 2)] + ["X"]:
            ccolors.append("#252525")
            chrs.add(c)

    D = 1 * (I_mtr.values == 1) * (O_mtr.values == 0) + 2 * (I_mtr.values == 0) * (
        O_mtr.values == 1
    )
    rvb = mcolors.LinearSegmentedColormap.from_list(
        "", ["#DEDEDE", "#D92347", "#A9D0F5", "#FFFFFF"]
    )

    # x = I.shape[0] / min(I.shape[0], I.shape[1])
    # y = I.shape[1] / min(I.shape[0], I.shape[1])

    return sns.clustermap(
        D,
        vmin=0,
        vmax=3,
        metric="euclidean",
        cmap=rvb,
        row_cluster=False,
        col_cluster=False,
        row_colors=[row_colors.values()],
        col_colors=[ccolors],
        cbar_pos=None,
        figsize=(11, 2),
        xticklabels=False,
        yticklabels=False,
        # linecolor='white',
        # linewidths=0.8,
        # annot=False,
        colors_ratio=(0.01, 0.05),
        dendrogram_ratio=0,
    )


def plot_noisy(df_in, row_colors=None):
    # TODO: check for test
    if row_colors is None:
        row_colors = {x: "#000000" for x in df_in.index}
    snvs = df_in.columns
    tmp = []
    for x in snvs:
        tmp.append(".".join(x.split(".chr")[1].split(".")[:2]))
    tmp = ns.index_natsorted(tmp)
    snvs = snvs[tmp]

    ccolors = []
    chrs = set()
    for x in snvs:
        c = x.split(".chr")[1].split(".")[0]
        if c in [f"{i}" for i in range(1, 22, 2)] + ["Y"]:
            ccolors.append("#969696")
            chrs.add(c)
        if c in [f"{i}" for i in range(2, 22, 2)] + ["X"]:
            ccolors.append("#252525")
            chrs.add(c)

    D = df_in.loc[row_colors.keys(), snvs]
    D.index.name = ""
    rvb = mcolors.LinearSegmentedColormap.from_list(
        "", ["#DEDEDE", "#000000", "#FFFFFF"]
    )

    return sns.clustermap(
        D,
        vmin=0,
        vmax=3,
        metric="euclidean",
        cmap=rvb,
        row_cluster=False,
        col_cluster=False,
        row_colors=[row_colors.values()],
        col_colors=[ccolors],
        cbar_pos=None,
        figsize=(8, 2),
        xticklabels=False,
        yticklabels=False,
        # linecolor='white',
        # linewidths=0.8,
        # annot=False,
        colors_ratio=(0.01, 0.05),
        dendrogram_ratio=0,
    )
