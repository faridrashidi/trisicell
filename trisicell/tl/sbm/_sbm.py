from collections import Counter

import numpy as np
import pandas as pd

import trisicell as tsc


def sbm(data2):
    graph_tool, graph_tool_is_not_imported = tsc.ul.import_graph_tool()
    if graph_tool_is_not_imported:
        raise RuntimeError("Unable to import a package!")

    from graph_tool.inference import minimize as gt_min

    def get_graph(df):
        G = graph_tool.Graph(directed=False)
        vtype = G.new_vertex_property("short")
        label2id = {}

        x = np.where(df == 1)
        for k in range(x[0].shape[0]):
            i = x[0][k]
            j = x[1][k]
            cell_key = f"cell{i}"
            mut_key = f"mut{j}"
            if cell_key in label2id:
                v = label2id[cell_key]
            else:
                v = G.add_vertex()
                label2id[cell_key] = int(v)
                vtype[v] = 1

            if mut_key in label2id:
                w = label2id[mut_key]
            else:
                w = G.add_vertex()
                label2id[mut_key] = int(w)
                vtype[w] = 2

            G.add_edge(v, w)
        return G, label2id, vtype

    def blockmodel_to_labels(b, label2id, n_cells=None, n_muts=None, maxblocks=100):
        if n_cells is None:
            n_cells = max(
                [int(x[4:]) if x.startswith("cell") else 0 for x in label2id.keys()]
            )
        if n_muts is None:
            n_muts = max(
                [int(x[3:]) if x.startswith("mut") else 0 for x in label2id.keys()]
            )

        cell_array = [
            (b[label2id[f"cell{i}"]] if f"cell{i}" in label2id else maxblocks + 1)
            for i in range(n_cells)
        ]
        temp = sorted(list(Counter(cell_array).keys()))
        cell_idx_to_blocknum = {temp[i]: i + 1 for i in range(len(temp))}
        cell_array = [cell_idx_to_blocknum[a] for a in cell_array]

        mut_array = [
            (b[label2id[f"mut{i}"]] if f"mut{i}" in label2id else maxblocks + 1)
            for i in range(n_muts)
        ]
        temp = sorted(list(Counter(mut_array).keys()))
        mut_idx_to_blocknum = {temp[i]: i + 1 for i in range(len(temp))}
        mut_array = [mut_idx_to_blocknum[a] for a in mut_array]
        return cell_array, mut_array

    def get_ratio(df):
        x = (df == 1).sum().sum() / df.size
        return x, 1 - x

    data = data2.copy()
    data[data == 3] = 0
    G, label2id, vtype = get_graph(data)

    min_blocks = 2
    max_blocks = 10
    r = gt_min.minimize_nested_blockmodel_dl(
        G, B_min=min_blocks, B_max=max_blocks, state_args={"clabel": vtype}
    )
    b = r.get_bs()[0]

    cell_array, mut_array = blockmodel_to_labels(
        b, label2id, n_cells=data.shape[0], n_muts=data.shape[1]
    )

    ixgrid = np.ix_(np.argsort(cell_array), np.argsort(mut_array))
    N = pd.DataFrame(
        data.values[ixgrid], index=np.sort(cell_array), columns=np.sort(mut_array)
    )
    I = data2.loc[
        data2.index[np.argsort(cell_array)], data2.columns[np.argsort(mut_array)]
    ]

    for n in set(N.columns):
        for m in set(N.columns):
            if get_ratio(N.loc[n, m])[0] > 0:
                N.loc[n, m] = 1

    N.index = I.index
    N.columns = I.columns
    out = N.loc[data2.index, data2.columns]

    return out
