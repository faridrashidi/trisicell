from collections import OrderedDict

import click

import trisicell as tsc
from trisicell.tl import consensus

from ._bnb import bnb
from ._booster import booster
from ._consensus import consensus
from ._huntress import huntress
from ._mcalling import mcalling
from ._phiscs import phiscsb, phiscsi
from ._scistree import scistree
from ._scite import scite
from ._score import score
from ._search import search
from ._trees import cf2newick, cf2tree


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.
    Make sure you initialize the `self.commands` with OrderedDict instance.
    With decorator, use::
        @click.group(cls=NaturalOrderGroup, commands=OrderedDict())
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.
        If the dict is OrderedDict, it will preserve the order commands
        were added.
        """
        return self.commands.keys()


@click.version_option(version=tsc.__version__)
@click.group(
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings=dict(max_content_width=300, terminal_width=300),
)
def cli():
    f"""Scalable intratumor heterogeneity inference and validation from single-cell data."""
    return None


cli.add_command(mcalling)
cli.add_command(score)
cli.add_command(scistree)
cli.add_command(scite)
cli.add_command(booster)
cli.add_command(phiscsb)
cli.add_command(phiscsi)
cli.add_command(bnb)
cli.add_command(huntress)
cli.add_command(cf2newick)
cli.add_command(cf2tree)
cli.add_command(consensus)
cli.add_command(search)
