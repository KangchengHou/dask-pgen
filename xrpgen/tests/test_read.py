import xrpgen
import numpy as np
import pandas_plink
import time


def test_phase():
    path = xrpgen.get_test_data("phased.pgen")
    d1 = xrpgen.read_pgen(path, phase=False)
    d2 = xrpgen.read_pgen(path, phase=True)
    assert np.all(d1 == d2.sum(axis=2)).compute()


def test_pandas_plink_consistency():
    path = xrpgen.get_test_data("plink1.bed")
    d1 = pandas_plink.read_plink1_bin(path, ref="a0", verbose=False).data
    d2 = xrpgen.read_pgen(path, phase=False)
    assert np.allclose(d1, d2, equal_nan=True)


def test_benchmark_pandas_plink(benchmark):
    path = xrpgen.get_test_data("plink1.bed")
    dset = pandas_plink.read_plink1_bin(path, ref="a0", verbose=False).data
    benchmark(dset.compute)


def test_benchmark_xrpgen(benchmark):
    path = xrpgen.get_test_data("plink1.bed")
    dset = xrpgen.read_pgen(path, phase=False)
    benchmark(dset.compute)