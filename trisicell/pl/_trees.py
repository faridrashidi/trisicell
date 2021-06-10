import networkx as nx
import numpy as np
import pandas as pd
from IPython.display import SVG, Image, display

import trisicell as tsc
from trisicell.pl._annotation import _add_barplot, _add_chromplot, _get_tree
from trisicell.ul._trees import _to_newick, to_mtree


def clonal_tree(
    tree,
    muts_as_number=False,
    cells_as_number=False,
    show_id=False,
    cell_info=None,
    output_file=None,
    color_attr=None,
    dpi=150,
):
    """Draw the tree in clonal format.

    This functions plots the tree in which edges are mutations and nodes are cells.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The input tree.
    muts_as_number : :obj:`bool`, optional
        Change the mutation list to a number at edges, by default False
    show_cells : :obj:`bool`, optional
        Whether to show cells or not, by default True
    collapsed_path : :obj:`bool`, optional
        Collapsed linear path in the tree , by default False
    cell_info : :class:`pandas.DataFrame`, optional
        Information of cells for coloring the nodes by a pie chart, by default None
    output_file : :obj:`str`, optional
        Path to a file for saving the tree in, by default None
    color_attr : :obj:`str`, optional
        Attributes in the `cell_info` dataframe for coloring the nodes, by default None
    dpi : :obj:`int`, optional
        Resolution of rendered figures – this influences the size of
        figures in notebooks, by default 150

    Returns
    -------
    :obj:`None`
    """

    graphviz, graphviz_is_not_imporeted = tsc.ul.import_graphviz()
    if graphviz_is_not_imporeted:
        raise RuntimeError("Unable to import a package!")

    tc = tree.copy()
    root = tsc.ul.root_id(tree)
    tc.nodes[root]["label"] = "root"
    tc.nodes[root]["fontname"] = "Helvetica"
    tc.nodes[root]["style"] = "rounded"
    tc.nodes[root]["shape"] = "box"
    tc.nodes[root]["margin"] = 0.05
    tc.nodes[root]["pad"] = 0
    tc.nodes[root]["width"] = 0
    tc.nodes[root]["height"] = 0

    if muts_as_number:
        for u, v, l in tc.edges.data("label"):
            ll = l.split(tc.graph["splitter_mut"])
            tc.add_edge(u, v, label=f"  {len(ll)}  ")

    if cells_as_number:
        for n in tc.nodes:
            if n != root:
                ll = tc.nodes[n]["label"].split(tc.graph["splitter_cell"])
                tc.nodes[n]["label"] = f"{len(ll)}"

    if cell_info is not None:
        tc.nodes[root]["label"] = tree.graph["splitter_cell"].join(
            tc.graph["normal_cells"]
        )
        mapping = cell_info[color_attr].to_dict()
        for node in tc:
            num = 0
            paths = nx.shortest_path(tc, source=root, target=node)
            for i in range(len(paths) - 1):
                x = paths[i]
                y = paths[i + 1]
                num += len(tc[x][y]["label"].split(tc.graph["splitter_mut"]))
            try:
                freq = [
                    mapping[x]
                    for x in tc.nodes[node]["label"].split(tc.graph["splitter_cell"])
                ]
            except:
                freq = ["#FFFFFF"]
            freq = pd.DataFrame(freq)[0].value_counts(normalize=True)
            fillcolor = ""
            for index, value in freq.items():
                fillcolor += f"{index};{value}:"
            tc.nodes[node]["fontsize"] = 14
            tc.nodes[node]["shape"] = "circle"
            tc.nodes[node]["fontname"] = "Helvetica"
            tc.nodes[node]["style"] = "wedged"
            tc.nodes[node]["margin"] = 0.05
            tc.nodes[node]["pad"] = 0
            if "––" in tc.nodes[node]["label"]:
                tc.nodes[node]["width"] = 0
                tc.nodes[node]["height"] = 0
            else:
                tc.nodes[node]["width"] = 0.8
                tc.nodes[node]["height"] = 0.8
            tc.nodes[node]["fillcolor"] = fillcolor  # "red;0.3:green;0.6:orange"
            tc.nodes[node]["label"] = len(
                tc.nodes[node]["label"].split(tc.graph["splitter_cell"])
            )
            tc.nodes[node]["fontcolor"] = "white"
            tc.nodes[node]["color"] = "gray"

    if show_id:
        for u, v, l in tc.edges.data("label"):
            tc.add_edge(u, v, label=l + f"\n[{v}]")
            tc.nodes[v]["label"] = tc.nodes[v]["label"] + f"\n[{v}]"

    tc.graph["graph"] = {"fontname": "Helvetica"}
    tc.graph["node"] = {"fontname": "Helvetica", "fontsize": 14}
    tc.graph["edge"] = {"fontname": "Helvetica", "fontsize": 14}

    tree.graph["type"] = "clonal"
    cell_list, mutation_list = _clonal_cell_mutation_list(tree)
    tree.graph["cell_list"] = cell_list
    tree.graph["mutation_list"] = mutation_list

    if output_file is not None:
        tsc.io.to_png(tc, output_file, dpi)

    return display(
        Image(nx.drawing.nx_pydot.to_pydot(tc).create_png(), embed=True, retina=True)
    )


def dendro_tree(
    tree,
    width=1200,
    height=500,
    dpi=300,
    cell_info=None,
    label_color="group_color",
    line_size=0.4,
    tiplab_size=2.5,
    inner_node_type="nmuts",
    inner_node_size=2,
    distance_labels_to_bottom=4,
    annotation=[],
    output_file=None,
):
    """Draw the tree in dendro fromat.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The input tree.
    width : :obj:`int`, optional
        Width of the figure, by default 8
    height : :obj:`int`, optional
        Height of the figure, by default 3
    cell_info : :class:`pandas.DataFrame`, optional
        Information about cells such as color,
        expression values of genes and etc, by default None
    label_color : :obj:`str`, optional
        The column name in which colors of cells are stored
        in the dataframe provided as `cell_info`, by default "group_color"
    line_size : :obj:`float`, optional
        Line size of the tree, by default 0.4
    tiplab_size : :obj:`float`, optional
        Cell name size in the tree, by default 2.5
    inner_node_type : :obj:`str`, optional
        The format of the inner nodes (i.e. mutations), by default "nmuts"
        Values are:

            - `nmuts`: only number of mutations
            - `nodeid`: only node id
            - `both`: both number of mutations and node id
    inner_node_size : :obj:`int`, optional
        Size of the inner nodes (i.e. mutations), by default 2
    distance_labels_to_bottom : :obj:`int`, optional
        Distance cell names to the bottom of the figure, by default 4
    annotation : :obj:`list`, optional
        List of gene names provided in the column dataframe of ``cell_info``
        in to be annotated in the bottom of the tree, by default []
    output_file : :obj:`str`, optional
        Path to a file for saving the tree in, by default None

    Returns
    -------
    :obj:`None`

    Note
    ----
    The cell names in the tree must be identical to the index of `cell_info`
    dataframe if it was provided.
    """

    if inner_node_type.lower() not in ["nmuts", "nodeid", "both"]:
        raise ValueError("Wrong `inner_node_type` choice!")

    ggtree, ggtree_is_not_imported = tsc.ul.import_rpy2(
        "ggtree",
        "devtools::install_github(c('YuLab-SMU/ggtree','xiangpin/ggtreeExtra','YuLab-SMU/aplot'))\n"
        "install.packages('cowplot')\n",
    )
    if ggtree_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.lib import grdevices
    from rpy2.robjects.packages import importr

    newick, info2, mutation_list = _newick_info2_mutation_list(tree)
    tree.graph["mutation_list"] = mutation_list.set_index("index")
    tree.graph["newick"] = newick
    tree.graph["type"] = "dendro"

    ggplot2 = importr("ggplot2")
    cowplot = importr("cowplot")
    ggtree = importr("ggtree")
    ape = importr("ape")
    aplot = importr("aplot")
    # ggtext = importr("ggtext")

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        if cell_info is not None:
            if "group_color" not in cell_info.columns:
                cell_info["group_color"] = tree.graph["data"].shape[0] * ["#000000"]
            cell_info_r = ro.conversion.py2rpy(cell_info.reset_index())
        else:
            cell_info_r = ro.conversion.py2rpy(
                pd.DataFrame.from_dict(
                    {
                        "x": tree.graph["data"].index,
                        "group_color": tree.graph["data"].shape[0] * ["#000000"],
                    }
                )
            )
        ro.globalenv["info1"] = cell_info_r
        info2_r = ro.conversion.py2rpy(info2)
        ro.globalenv["info2"] = info2_r

    # p = ggtree.ggtree(
    #     ape.read_tree(text=f"{newick}"), layout="dendrogram", size=line_size
    # )
    cmd = _get_tree(
        newick,
        line_size,
        label_color,
        tiplab_size,
        inner_node_type,
        inner_node_size,
        distance_labels_to_bottom,
    )

    for i, ann in enumerate(annotation):
        if i == 0:
            cmd += "info1$id_index <- p$data[match(info1[,1], p$data$label),]$y"
        if ann[0] == "bar":
            cmd += _add_barplot(ann[1], ann[2], ann[3])
        elif ann == "chrom":
            cmd += _add_chromplot()
        elif ann == "imputation":
            for imp in [
                ("f_0_1", "f_0_1_color"),
                ("f_1_0", "f_1_0_color"),
                ("f_3_0", "f_3_0_color"),
                ("f_3_1", "f_3_1_color"),
            ]:
                cmd += _add_barplot(imp[0], imp[1])

    with grdevices.render_to_bytesio(
        grdevices.png, width=width, height=height, res=dpi
    ) as image:
        p = ro.r(cmd)
        ro.r.show(p)
        if output_file is not None:
            ro.r.ggsave(
                plot=p,
                filename=output_file,
                width=width / dpi,
                height=height / dpi,
                units="in",
                dpi=dpi,
                limitsize=False,
            )
    return display(Image(image.getvalue(), embed=True, retina=True))


def _clonal_cell_mutation_list(tree):
    muts_list = []
    cells_list = []
    for u, v, l in tree.edges.data("label"):
        muts = l.split(tree.graph["splitter_mut"])
        if "––" not in tree.nodes[v]["label"]:
            cells = tree.nodes[v]["label"].split(tree.graph["splitter_cell"])
            for mut in muts:
                muts_list.append({"mut": mut, "Node": f"[{v}]"})
            for cell in cells:
                cells_list.append({"cell": cell, "Node": f"[{v}]"})
        else:
            for mut in muts:
                muts_list.append({"mut": mut, "Node": f"[{v}]"})
    cells_list = pd.DataFrame(cells_list).set_index("cell")
    muts_list = pd.DataFrame(muts_list).set_index("mut")
    return cells_list, muts_list


def _newick_info2_mutation_list(tree):
    tree2 = to_mtree(tree)

    row = []
    for node in tree2.nodes:
        if tree2.in_degree(node) == 0:
            row.append(
                {
                    "newick_label": f"Node{node+1}",
                    "nmuts_label": "root",
                    "nodeid_label": f"[{node+1}]",
                    "both_label": f"[{node+1}]: root",
                }
            )
        else:
            row.append(
                {
                    "newick_label": f"Node{node+1}",
                    "nmuts_label": f"{len(tree2.nodes[node]['label'])}",
                    "nodeid_label": f"[{node+1}]",
                    "both_label": f"[{node+1}]: {len(tree2.nodes[node]['label'])}",
                }
            )
    info2 = pd.DataFrame(row)

    row = []
    for node in tree2.nodes:
        if tree2.in_degree(node) == 0:
            continue
        for mut in tree2.nodes[node]["label"]:
            ens, gene, chrom, pos, ref, alt = tsc.ul.split_mut(mut)
            if ens is None:
                row.append(
                    {
                        "index": mut,
                        "node_id": f"[{node+1}]",
                    }
                )
            else:
                row.append(
                    {
                        "index": mut,
                        "Node": f"[{node+1}]",
                        "Ensemble": ens,
                        "Gene": gene,
                        "Chrom": chrom,
                        "Position": pos,
                        "Reference": ref,
                        "Alteration": alt,
                    }
                )
    mutation_list = pd.DataFrame(row)
    newick = _to_newick(tree)
    return newick, info2, mutation_list
