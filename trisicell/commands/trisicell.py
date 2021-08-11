from collections import OrderedDict

import click

import trisicell as tsc
from trisicell.commands._bnb import bnb
from trisicell.commands._booster import booster
from trisicell.commands._consensus import consensus
from trisicell.commands._huntress import huntress
from trisicell.commands._mcalling import mcalling
from trisicell.commands._phiscs import phiscsb, phiscsi
from trisicell.commands._scistree import scistree
from trisicell.commands._scite import scite
from trisicell.commands._score import score
from trisicell.commands._search import search
from trisicell.commands._trees import cf2newick, cf2tree


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


# def citation(ctx, param, value):
#     print("Please check https://trisicell.readthedocs.io/citing.html")
#     ctx.exit(0)


@click.version_option(version=tsc.__version__)
# @click.option("--cite", is_flag=True, callback=citation, help="Show citation bib.")
@click.group(
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings={"max_content_width": 300, "terminal_width": 300},
)
def cli():
    """Trisicell.

    Scalable intratumor heterogeneity inference and validation from single-cell data.
    """
    return None


cli.add_command(mcalling)
cli.add_command(booster)
cli.add_command(phiscsb)
cli.add_command(phiscsi)
cli.add_command(scite)
cli.add_command(scistree)
cli.add_command(bnb)
cli.add_command(huntress)
cli.add_command(cf2newick)
cli.add_command(cf2tree)
cli.add_command(score)
cli.add_command(consensus)
cli.add_command(search)
