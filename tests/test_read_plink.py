import dapgen
import numpy as np
import pandas_plink

# Test `dapgen.read_plink` function for multi-purpose use.
def test_basic():
    path = dapgen.get_test_data("plink2.merged.pgen")
    d = dapgen.read_plink(path, phase=False)


def test_consistent():
    path1 = dapgen.get_test_data("plink1.merged.bed")
    path2 = dapgen.get_test_data("plink2.merged.pgen")
    d1 = pandas_plink.read_plink1_bin(path1, ref="a0", verbose=False).data.T
    geno1, df_snp1, df_indiv1 = dapgen.read_plink(path1, phase=False)
    geno2, df_snp2, df_indiv2 = dapgen.read_plink(path2, phase=False)
    assert np.allclose(geno1, geno2)
    assert df_snp1.index.equals(df_snp2.index)
    assert np.all(df_indiv1.IID.values == df_indiv2.index.values)


def test_two_chrom():
    path1 = dapgen.get_test_data("plink2.chr21.pgen")
    path2 = dapgen.get_test_data("plink2.chr22.pgen")
    merged_path = dapgen.get_test_data("plink2.merged.pgen")
    geno1, df_snp1, df_indiv1 = dapgen.read_plink(merged_path, phase=True)

    for p in [[path1, path2], [path2, path1]]:
        # read_plink can also infer the order of chromosomes
        geno2, df_snp2, df_indiv2 = dapgen.read_plink(p, phase=True)
        assert np.allclose(geno1, geno2)
        assert df_snp1.equals(df_snp2)
        assert df_indiv1.equals(df_indiv2)


def test_two_indiv():
    path1 = dapgen.get_test_data("plink2.merged.FIN.pgen")
    path2 = dapgen.get_test_data("plink2.merged.GBR.pgen")
    merged_path = dapgen.get_test_data("plink2.merged.pgen")
    geno1, df_snp1, df_indiv1 = dapgen.read_plink(merged_path, phase=True)

    # read_plink can also infer the order of chromosomes
    geno2, df_snp2, df_indiv2 = dapgen.read_plink([path1, path2], phase=True)
    assert df_snp1.equals(df_snp2)
    for pop in ["GBR", "FIN"]:
        assert np.allclose(
            geno1[:, df_indiv1.Population == pop],
            geno2[:, df_indiv2.Population == pop],
        )