import matplotlib.colors as mcolors
import seaborn as sns


def heatmap(adata, color_attrs=None, layer="X", figsize=(12, 7)):
    # TODO: check for test
    if color_attrs is not None:
        row_colors = []
        for attr in color_attrs:
            row_colors.append(adata.obs[attr])
    else:
        row_colors = None

    # if c in [f"{i}" for i in range(1, 22, 2)] + ["Y"]:
    #     ccolors.append("#969696")
    #     chrs.add(c)
    # if c in [f"{i}" for i in range(2, 22, 2)] + ["X"]:
    #     ccolors.append("#252525")
    #     chrs.add(c)

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
