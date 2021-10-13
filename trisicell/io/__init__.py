"""Input/Output Module."""

from trisicell.io._genotype import read, write
from trisicell.io._tree import to_png

__all__ = (
    read,
    write,
    to_png,
)
