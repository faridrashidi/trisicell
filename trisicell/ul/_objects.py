import pandas as pd
from IPython import get_ipython


class Cell:
    def __init__(self, n_chroms=3):
        super().__init__()
        self.n_chroms = n_chroms
        self.image = "/data/frashidi/test.png"

    def __str__(self):
        """Print."""
        return f"n_chroms={self.n_chroms}"

    def plot(self, path):
        ipython = get_ipython()
        ipython.magic("matplotlib agg")
        df = pd.DataFrame({"lab": ["A", "B", "C"], "val": [10, 30, 20]})
        ax = df.plot.bar(x="lab", y="val", rot=0)
        ax.get_figure().savefig(path, bbox_inches="tight", pad_inches=0)
        self.image = path
        ipython.magic("matplotlib inline")
