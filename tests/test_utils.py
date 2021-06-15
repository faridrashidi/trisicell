import numpy as np
import pandas as pd

import trisicell as tsc


class TestUtils:
    def test_hclustering(self):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randint(0, 2, size=(10, 10)))
        clusters = tsc.ul.hclustering(df)
        assert clusters[6].value_counts()[5] == 3
