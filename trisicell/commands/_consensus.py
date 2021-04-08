import click

import trisicell as tsc


@click.command(short_help="Calculate consensus betweeen two trees.")
@click.argument(
    "first_tree",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "second_tree",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "consensus_path",
    required=True,
    type=click.Path(exists=False),
)
def consensus(first_tree, second_tree, consensus_path):
    """The consensus tree between two phylogenetic trees.

    It writes the conflict-free matrix representing the consensus tree
    into the `consensus_path` filepath.

    trisicell consensus first_tree.CFMatrix second_tree.CFMatrix
    """

    tsc.settings.verbosity = "info"

    sc1 = tsc.io.read(first_tree)
    sc2 = tsc.io.read(second_tree)
    final_tree = tsc.tl.consensus_tree(sc1, sc2)
    tsc.io.write(final_tree.graph["data"], consensus_path)

    return None
