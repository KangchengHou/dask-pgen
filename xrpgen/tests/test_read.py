import xrpgen
import numpy as np


def test_phase():
    d1 = xrpgen.read_pgen(xrpgen.get_test_data("phased.pgen"), phase=False)
    d2 = xrpgen.read_pgen(xrpgen.get_test_data("phased.pgen"), phase=True)
    assert np.all(d1 == d2.sum(axis=2)).compute()