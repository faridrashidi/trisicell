import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import natsort as ns
import numpy as np
import pandas as pd
import seaborn as sns

import trisicell as tsc


def heatmap(adata, color_attrs, layer="X", figsize=(12, 7)):
    row_colors = []
    for attr in color_attrs:
        row_colors.append(adata.obs[attr])

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
        # vmin=0,
        # vmax=1,
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


def draw(adata, figsize=(10, 2), dpi=100):
    gspec = plt.GridSpec(1, 2, plt.figure(None, figsize, dpi=dpi))
    ax1 = plt.subplot(gspec[0])

    ax1 = plt.subplot(gspec[1])
    ax2.margins(0)

    return [ax1, ax2]


def plot_size(df_in):
    plt.rcParams["xtick.labelsize"] = 6
    plt.rcParams["ytick.labelsize"] = 6
    plt.figure(figsize=(1, 1))
    ax = sns.histplot([], color="black")
    ax.set_xlabel(f"{df_in.shape[1]} mutations", fontsize=8)
    ax.set_ylabel(f"{df_in.shape[0]} cells", fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.xaxis.set_label_position("top")


def plot_input(df_in):
    plt.rcParams["xtick.labelsize"] = 6
    plt.rcParams["ytick.labelsize"] = 6
    plt.figure(figsize=(1, 1))
    df = pd.DataFrame(np.array(np.unique(df_in.values, return_counts=True)).T)
    df[2] = 100 * df[1] / df_in.size
    df = df.set_index([0]).rename({0: "0", 1: "1", 3: "NA"})
    if "NA" in df.index:
        ax = df[2].plot.barh(color="#D92347", width=0.7)
    else:
        ax = df[2].plot.barh(color="#D92347", width=0.5)
    plt.gca().invert_yaxis()
    plt.xlim(0, 100)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks(np.arange(0, 101, 25))
    ax.tick_params(axis="x", direction="out", length=2, width=1)
    ax.tick_params(axis="y", direction="out", length=2, width=1)


def plot_output(df_out):
    plt.rcParams["xtick.labelsize"] = 6
    plt.rcParams["ytick.labelsize"] = 6
    plt.figure(figsize=(1, 1))
    df = pd.DataFrame(np.array(np.unique(df_out.values, return_counts=True)).T)
    df[2] = 100 * df[1] / df_out.size
    df = df.set_index([0]).rename({0: "0", 1: "1", 3: "NA"})
    ax = df[2].plot.barh(color="#0067A5", width=0.5)
    plt.gca().invert_yaxis()
    plt.xlim(0, 100)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks(np.arange(0, 101, 25))
    ax.tick_params(axis="x", direction="out", length=2, width=1)
    ax.tick_params(axis="y", direction="out", length=2, width=1)


def plot_flips2(df_in, df_out):
    I = df_in
    O = df_out
    df = pd.DataFrame(0, index=["0->1", "1->0", "NA->0", "NA->1"], columns=[1])
    df.loc["0->1", 1] = ((I == 0) & (O == 1)).sum().sum()
    df.loc["1->0", 1] = ((I == 1) & (O == 0)).sum().sum()
    df.loc["NA->0", 1] = ((I == 3) & (O == 0)).sum().sum()
    df.loc["NA->1", 1] = ((I == 3) & (O == 1)).sum().sum()
    print(df)


def plot_flips(df_in, df_out, row_colors):
    snvs = np.intersect1d(df_in.columns, df_out.columns)
    # cells = np.intersect1d(df_in.index, df_out.index)
    tmp = []
    for x in snvs:
        tmp.append(".".join(x.split(".chr")[1].split(".")[:2]))
    tmp = ns.index_natsorted(tmp)
    snvs = snvs[tmp]
    I = df_in.loc[row_colors.keys(), snvs]
    O = df_out.loc[row_colors.keys(), snvs]

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

    D = 1 * (I.values == 1) * (O.values == 0) + 2 * (I.values == 0) * (O.values == 1)
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
