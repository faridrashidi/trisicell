import sys

__author__ = "Simone Ciccolella"
__date__ = "11/30/21"

name = "nh2lgf"

desc = """

convert newhampshire (or newick) format (nh)
to lemon graph format (lgf)
-- Murray Patterson, Sept 2014

"""

note = ""


# for processing a 'node' object, i,e., of the form 'A:B', where 'A'
# is a name (or a bootstrap value in the case of an internal node),
# and 'B' is a branch length
def node(suffix, count, nodes, arcs):

    assert not suffix.startswith("("), "\n\ninput not in nh format"

    # eat until ',', ')' or ';'
    if ")" not in suffix:  # we know we're at the end
        a, b, c = suffix.partition(";")  # if not, search for ','
    elif "," in suffix and suffix.find(",") < suffix.find(")"):
        a, b, c = suffix.partition(",")
    else:  # o.w., it is the only terminal
        a, b, c = suffix.partition(")")

    name, x, brlen = a.partition(":")
    nodes.append((count + 1, name))

    return b + c, count + 1, brlen, nodes, arcs


# for processing a 'subtree' object, i.e., of the form 'N' or 'X',
# where 'N' is a 'node' object (base case), and 'X' is a 'subtrees'
# object (recursive case)
def subtree(suffix, count, nodes, arcs):

    if suffix.startswith("("):  # recursive case
        return subtrees(suffix, count, nodes, arcs)

    else:  # base case, we're at a node (which is necessarily a leaf, here)
        return node(suffix, count, nodes, arcs)


# for recieving and processing a 'subtrees' object, i.e. of the form
# '(S,S,..,S)' where 'S' is a 'subtree' object
def subtrees(suffix, count, nodes, arcs):

    assert ")" in suffix, "\n\ninput not in nh format"

    suffix = suffix[1:]  # eat the '('
    children = []  # store info for each 'S' in '(S,S,..,S)'

    # loop through the 'S' in '(S,S,..,S)'
    while True:  # note : contains at least one 'S', even if ''
        suffix, count, brlen, nodes, arcs = subtree(suffix, count, nodes, arcs)
        children.append((count, brlen))

        if suffix.startswith(")"):  # stopping condition
            break
        if suffix.startswith(","):  # in case of multiple subtrees
            suffix = suffix[1:]  # eat the ','

    # now, here, suffix[0] == ')'
    suffix = suffix[1:]  # eat the ')'
    # get parent of '(S,S,..,S)'
    suffix, count, brlen, nodes, arcs = node(suffix, count, nodes, arcs)
    for child in children:  # now add the arcs to each child 'S'
        arcs.append((child[0], count, child[1]))

    return suffix, count, brlen, nodes, arcs


def newick_to_edgelist(input_string):
    s = input_string
    assert s.find(";") == len(s) - 1, "\n\nerror: input not in nh format"
    assert s.count("(") == s.count(")"), "\n\nerror: input not in nh format"

    # parse the input
    nodes = []
    arcs = []
    subtree(s, -1, nodes, arcs)

    try:
        # output in lgf format
        node_dict = {}
        for label, name in nodes:
            node_dict[label] = name

        edges = []
        for u, v, _ in arcs:
            edges.append((u, v))

        return node_dict, edges
    except OSError:
        sys.exit(0)
