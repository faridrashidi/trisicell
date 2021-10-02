import networkx as nx
from IPython.display import SVG, Image, display

import trisicell as tsc


def simulate(tree, show_image=False, kind="png"):
    tmpdir = tsc.ul.tmpdirsys(suffix=".simulate")
    tree2 = tree.copy()
    for n in tree2.nodes():
        node = tree2.nodes[n]
        if show_image:
            node["cell"].plot(f"{tmpdir.name}/salam.png")
            node["image"] = node["cell"].image
            node["shape"] = None
            node["color"] = None
            node["margin"] = 0
            node["width"] = 0
            node["height"] = 0
    if kind == "png":
        display(Image(nx.drawing.nx_pydot.to_pydot(tree2).create_png(), retina=True))
    elif kind == "svg":
        display(SVG(nx.drawing.nx_pydot.to_pydot(tree2).create_svg()))
    tmpdir.cleanup()
