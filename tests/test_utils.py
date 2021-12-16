import pandas as pd
import dapgen
import numpy as np
from typing import List


def test_align_snp():
    df1 = pd.DataFrame(
        {
            "CHROM": [1, 1, 1, 1, 1],
            "POS": [1, 2, 3, 4, 5],
            "REF": ["A", "A", "A", "A", "A"],
            "ALT": ["T", "T", "T", "T", "T"],
        },
        index=[0, 1, 2, 3, 4],
    )
    df2 = pd.DataFrame(
        {
            "CHROM": [1, 1, 1, 1],
            "POS": [1, 2, 4, 5],
            "REF": ["A", "T", "A", "A"],
            "ALT": ["T", "A", "T", "T"],
        },
        index=[0, 1, 2, 3],
    )
    idx1, idx2, flip_sign = dapgen.align_snp(df1, df2)
    assert np.allclose(idx1.values, [0, 1, 3, 4])
    assert np.allclose(idx2.values, [0, 1, 2, 3])
    assert np.allclose(flip_sign, [1.0, -1.0, 1.0, 1.0])
