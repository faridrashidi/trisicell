import copy
import os

import networkx as nx
import numpy as np
import pandas as pd
from IPython.display import SVG, Image, display

import trisicell as tsc
from trisicell.pl._annotation import _add_barplot, _add_chromplot, _get_tree
from trisicell.ul._trees import _info2_mutation_list, _to_newick


def clonal_tree(
    tree,
    muts_as_number=False,
    show_cells=True,
    collapsed_path=False,
    cell_info=None,
    output_file=None,
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

    Returns
    -------
    :obj:`None`
    """

    tc = tree.copy()
    root = None
    for node in tc.nodes:
        if tc.in_degree(node) == 0:
            root = node
            break

    if collapsed_path:
        muts_as_number = True
        show_cells = True

    if muts_as_number:
        for u, v, l in tc.edges.data("label"):
            ll = l.split(tc.graph["splitter_mut"])
            if "––" in tc.nodes[v]["label"]:
                tc.add_edge(
                    u, v, label=f"  {len(ll)}  ", color="black", fontcolor="black"
                )
            else:
                tc.add_edge(u, v, label=f"  {len(ll)}  ")

    if not show_cells:
        for node in tc.nodes:
            if node != root:
                tc.nodes[node]["label"] = ""
                tc.nodes[node]["shape"] = "circle"
                tc.nodes[node]["width"] = 0
                tc.nodes[node]["height"] = 0

    if collapsed_path:
        tc2 = tc.copy()
        for _ in range(len(tc.nodes)):
            d_in = tc2.in_degree(tc2)
            d_out = tc2.out_degree(tc2)
            for node in tc2.nodes():
                if d_out[node] == 1 and d_in[node] == 1:
                    parent = [x for x in tc2.predecessors(node)][0]
                    child = [x for x in tc2.successors(node)][0]
                    if d_out[parent] < 2 and d_in[parent] == 1:
                        new_node = f"{parent}+{node}"
                        new_label = (
                            f"{tc2.nodes[parent]['label']}"
                            f"{tc.graph['splitter_cell']}"
                            f"{tc2.nodes[node]['label']}"
                        )
                        a = int(tc2[parent][node]["label"])
                        b = int(tc2[node][child]["label"])
                        new_edge = f"{a + b}"

                        tc2 = nx.contracted_nodes(tc2, parent, node, self_loops=False)
                        mapping = {parent: new_node}
                        tc2 = nx.relabel_nodes(tc2, mapping)
                        tc2[new_node][child]["label"] = new_edge
                        tc2.nodes[new_node]["label"] = new_label
                        break
        d_in = tc2.in_degree(tc2)
        d_out = tc2.out_degree(tc2)
        nodes = []
        for node in tc2.nodes():
            if d_out[node] == 0:
                nodes.append(node)
        for node in nodes:
            parent = [x for x in tc2.predecessors(node)][0]
            if d_out[parent] == 1 and d_in[parent] == 1:
                grandparent = [x for x in tc2.predecessors(parent)][0]

                new_node = f"{parent}+{node}"
                new_label = (
                    f"{tc2.nodes[parent]['label']}"
                    f"{tc.graph['splitter_cell']}"
                    f"{tc2.nodes[node]['label']}"
                )
                a = int(tc2[grandparent][parent]["label"])
                b = int(tc2[parent][node]["label"])
                new_edge = f"{a + b}"

                tc2 = nx.contracted_nodes(tc2, parent, node, self_loops=False)
                mapping = {parent: new_node}
                tc2 = nx.relabel_nodes(tc2, mapping)
                tc2[grandparent][new_node]["label"] = new_edge
                tc2.nodes[new_node]["label"] = new_label
        tc = tc2

        mapping = cell_info["name"].to_dict()
        for node in tc:
            if node != root:
                num = 0
                paths = nx.shortest_path(tc, source=root, target=node)
                for i in range(len(paths) - 1):
                    x = paths[i]
                    y = paths[i + 1]
                    num += len(tc[x][y]["label"].split(tc.graph["splitter_mut"]))
                try:
                    freq = [
                        mapping[x].split("_")[1]
                        for x in tc.nodes[node]["label"].split(
                            tc.graph["splitter_cell"]
                        )
                    ]
                except:
                    freq = ["#FFFFFF"]
                freq = pd.DataFrame(freq)[0].value_counts(normalize=True)
                fillcolor = ""
                for index, value in freq.items():
                    fillcolor += f"{index};{value}:"
                tc.nodes[node]["fontsize"] = 0.5
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
                tc.nodes[node]["label"] = node
            else:
                tc.nodes[node]["label"] = f"<<b>zygote</b>>"

    if output_file is not None:
        mygraph = nx.drawing.nx_agraph.to_agraph(tc)
        mygraph.layout(prog="dot")
        mygraph.draw(output_file)

    return display(SVG(nx.drawing.nx_pydot.to_pydot(tc).create_svg()))


def _merge(x):
    return (
        f"{x['Ensemble']}_{x['Gene']}.{x['Chrom']}."
        f"{x['Position']}.{x['Reference']}.{x['Alteration']}"
    )


def _get_newick_info2_mutations(tree):
    newick = _to_newick(tree)
    info2, mutation_list = _info2_mutation_list(tree)
    return newick, info2, mutation_list


def dendro_tree(
    tree,
    width=8,
    height=3,
    cell_info=None,
    label_color="group_color",
    line_size=0.4,
    tiplab_size=2.5,
    inner_node_id=False,
    inner_node_size=2,
    distance_labels_to_bottom=4,
    annotation=[],
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
    inner_node_id : :obj:`bool`, optional
        Information of the inner nodes (i.e. mutations), by default False
        Values are:

            - `nmuts_label`: only number of mutations
            - `nodeid_label`: only node id
            - `mixed_label`: both number of mutations and node id
    inner_node_size : :obj:`int`, optional
        Size of the inner nodes (i.e. mutations), by default 2
    distance_labels_to_bottom : :obj:`int`, optional
        Distance cell names to the bottom of the figure, by default 4
    annotation : :obj:`list`, optional
        List of gene names provided in the column dataframe of ``cell_info``
        in to be annotated in the bottom of the tree, by default []

    Returns
    -------
    :obj:`None`

    Raises
    ------
    SystemExit
        If `rpy2 <https://pypi.org/project/rpy2/>`_ is not already installed.
    SystemExit
        If `ggtree <https://github.com/YuLab-SMU/ggtree>`_ is not already installed.
    """

    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.lib import grdevices
        from rpy2.robjects.packages import importr
    except:
        raise SystemExit(
            "A problem was encountered importing `rpy2`. "
            "To run this `rpy2` and `R` need to be installed."
        )
    try:
        importr("ggtree")
    except:
        raise SystemExit(
            "A problem was encountered importing `ggtree` in R. "
            "To run this `ggtree` needs to be installed in R. "
            "Use the following lines to installed them.\n\n"
            "devtools::install_github(c('YuLab-SMU/ggtree',"
            "'xiangpin/ggtreeExtra','YuLab-SMU/aplot'))\n"
            "install.packages('cowplot')\n"
        )

    newick, info2, mutation_list = _get_newick_info2_mutations(tree)
    mutation_list["index"] = mutation_list.apply(_merge, axis=1)
    tree.graph["mutation_list"] = mutation_list.set_index("index")
    tree.graph["newick"] = newick

    ggplot2 = importr("ggplot2")
    cowplot = importr("cowplot")
    ggtree = importr("ggtree")
    ape = importr("ape")
    aplot = importr("aplot")
    pandas2ri.activate()

    if cell_info is not None:
        if "group_color" not in cell_info.columns:
            cell_info["group_color"] = tree.graph["data"].shape[0] * ["#000000"]
        cell_info_r = robjects.conversion.py2rpy(cell_info.reset_index())
    else:
        cell_info_r = robjects.conversion.py2rpy(
            pd.DataFrame.from_dict(
                {
                    "x": tree.graph["data"].index,
                    "group_color": tree.graph["data"].shape[0] * ["#000000"],
                }
            )
        )
    robjects.globalenv["info1"] = cell_info_r
    info2_r = robjects.conversion.py2rpy(info2)
    robjects.globalenv["info2"] = info2_r

    # p = ggtree.ggtree(
    #     ape.read_tree(text=f"{newick}"), layout="dendrogram", size=line_size
    # )
    cmd = _get_tree(
        newick,
        line_size,
        label_color,
        tiplab_size,
        inner_node_id,
        inner_node_size,
        distance_labels_to_bottom,
    )

    for ann in annotation:
        cmd += "info1$id_index <- p$data[match(info1[,1], p$data$label),]$y"
        if ann[0] == "bar":
            cmd += _add_barplot(ann[1], ann[2])
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
        grdevices.svg, width=width, height=height
    ) as image:
        p = robjects.r(cmd)
        robjects.r.show(p)
        # robjects.r.ggsave(
        #     plot=p,
        #     filename="x.pdf",
        #     device="pdf",
        #     width=7.5,
        #     height=3,
        #     units="in",
        #     limitsize=False,
        # )
    return display(SVG(image.getvalue()))
