import dapgen
import numpy as np
import pandas_plink


def test_phase():
    path = dapgen.get_test_data("plink2.merged.pgen")
    d1 = dapgen.read_pgen(path, phase=False)
    d2 = dapgen.read_pgen(path, phase=True)
    assert np.all(d1 == d2.sum(axis=2)).compute()


def test_pandas_plink_consistency():
    path = dapgen.get_test_data("plink1.merged.bed")
    d1 = pandas_plink.read_plink1_bin(path, ref="a0", verbose=False).data.T
    d2 = dapgen.read_pgen(path, phase=False)
    assert np.allclose(d1, d2)


def test_benchmark_pandas_plink(benchmark):
    path = dapgen.get_test_data("plink1.merged.bed")
    dset = pandas_plink.read_plink1_bin(path, ref="a0", verbose=False).data
    benchmark(dset.compute)


def test_benchmark_dapgen_plink1(benchmark):
    path = dapgen.get_test_data("plink1.merged.bed")
    dset = dapgen.read_pgen(path, phase=False)
    benchmark(dset.compute)


def test_benchmark_dapgen_plink2_nophase(benchmark):
    path = dapgen.get_test_data("plink2.merged.pgen")
    dset = dapgen.read_pgen(path, phase=False)
    benchmark(dset.compute)


def test_benchmark_dapgen_plink2_phase(benchmark):
    path = dapgen.get_test_data("plink2.merged.pgen")
    dset = dapgen.read_pgen(path, phase=True)
    benchmark(dset.compute)


# def test_benchmark_dapgen_plink2_no_phase(benchmark):
#     path = dapgen.get_test_data("1000G.EUR.QC.22.pgen")
#     dset = dapgen.read_pgen(path, phase=False)
#     benchmark(dset[:, 0:5000].compute)


# def test_benchmark_dapgen_plink2_phase(benchmark):
#     path = dapgen.get_test_data("1000G.EUR.QC.22.pgen")
#     dset = dapgen.read_pgen(path, phase=True)
#     benchmark(dset[:, 0:5000, :].compute)


# def test_benchmark_dapgen_plink2_no_phase(benchmark):
#     path = "/u/project/pasaniuc/pasaniucdata/admixture/projects/PAGE-QC/s01_vcf_dataset/hm3/merged.pgen"
#     dset = dapgen.read_pgen(path, phase=False)
#     benchmark(dset[0:5000, :].compute)


# def test_benchmark_dapgen_plink2_phase(benchmark):
#     path = "/u/project/pasaniuc/pasaniucdata/admixture/projects/PAGE-QC/s01_vcf_dataset/hm3/merged.pgen"
#     dset = dapgen.read_pgen(path, phase=True)
#     benchmark(dset[0:5000, :, :].compute)