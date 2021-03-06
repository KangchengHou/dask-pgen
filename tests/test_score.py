import pandas as pd
import dapgen
import numpy as np
from typing import List
import shutil
import os


def python_score(
    plink_path: str,
    df_weight: pd.DataFrame,
    weight_cols: List[str],
    center: bool = False,
    freq: np.ndarray = None,
):
    """
    Alternative way of implementing the score with python API, used for testing
    """
    geno, df_snp, df_indiv = dapgen.read_plink(plink_path)
    df_weight_info = ["CHROM", "POS", "REF", "ALT"]
    df_snp_idx, df_weight_idx, flip_sign = dapgen.align_snp(
        df_snp, df_weight[df_weight_info]
    )
    geno = geno[df_snp.index.get_indexer(df_snp_idx)].compute()
    df_snp = df_snp.loc[df_snp_idx]

    df_weight = df_weight.loc[df_weight_idx, weight_cols]
    to_flip = flip_sign == -1
    geno[to_flip, :] = 2.0 - geno[to_flip, :]
    if freq is None:
        center_values = geno.mean(axis=1)
    else:
        center_values = freq * 2
    if center:
        geno -= center_values[:, None]
    true_score = geno.T.dot(df_weight.values)
    return true_score, df_snp


def test_basic():
    # build df_weight
    plink_path = dapgen.get_test_data("plink2.merged.pgen")
    np.random.seed(1234)
    df_weight = dapgen.read_pvar(plink_path.replace(".pgen", ".pvar"))
    weight_cols = ["WEIGHT1", "WEIGHT2", "WEIGHT3"]
    for col in weight_cols:
        df_weight[col] = np.random.normal(size=len(df_weight))

    # full observation
    df_score, df_score_snp = dapgen.score(plink_path, df_weight)
    true_score, true_df_snp = python_score(plink_path, df_weight, weight_cols)
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)

    # test center
    df_score, df_score_snp = dapgen.score(plink_path, df_weight, center=True)
    true_score, true_df_snp = python_score(
        plink_path, df_weight, weight_cols, center=True
    )
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)

    # plink1
    plink_path = dapgen.get_test_data("plink1.merged.bed")
    df_score, df_score_snp = dapgen.score(plink_path, df_weight)
    true_score, true_df_snp = python_score(plink_path, df_weight, weight_cols)
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)

    # chrom 21 only
    plink_path = dapgen.get_test_data("plink2.chr21.pgen")
    df_score, df_score_snp = dapgen.score(plink_path, df_weight)
    true_score, true_df_snp = python_score(plink_path, df_weight, weight_cols)
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)

    # chrom 21 only with some snp flip in df_weight
    plink_path = dapgen.get_test_data("plink2.chr21.pgen")
    df_weight_flipped = df_weight.copy()
    flip_snp_idx = df_weight.index[[1, 3, 5]]
    for idx in flip_snp_idx:
        df_weight_flipped.loc[idx, ["REF", "ALT"]] = df_weight_flipped.loc[
            idx, ["ALT", "REF"]
        ].values
    df_score, df_score_snp = dapgen.score(plink_path, df_weight_flipped)
    true_score, true_df_snp = python_score(plink_path, df_weight_flipped, weight_cols)
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)

    # chrom 21 only with some wrongly coded allele
    plink_path = dapgen.get_test_data("plink2.chr21.pgen")
    df_weight_flipped = df_weight.copy()
    flip_snp_idx = df_weight.index[[1, 3, 5]]
    for idx in flip_snp_idx:
        df_weight_flipped.loc[idx, ["REF", "ALT"]] = df_weight_flipped.loc[
            idx, ["ALT", "REF"]
        ].values
    df_weight_flipped.loc[flip_snp_idx[0], "ALT"] = "WRONG"

    df_score, df_score_snp = dapgen.score(plink_path, df_weight_flipped)
    true_score, true_df_snp = python_score(plink_path, df_weight_flipped, weight_cols)
    assert np.allclose(df_score.values, true_score)
    assert df_score_snp.equals(true_df_snp)


def test_read_freq():
    # build df_weight
    plink_path = dapgen.get_test_data("plink2.merged.pgen")
    np.random.seed(1234)
    df_weight = dapgen.read_pvar(plink_path.replace(".pgen", ".pvar"))
    weight_cols = ["WEIGHT1", "WEIGHT2", "WEIGHT3"]
    for col in weight_cols:
        df_weight[col] = np.random.normal(size=len(df_weight))

    # full observation
    df_score, df_score_snp = dapgen.score(
        plink_path, df_weight, freq_suffix=".afreq", center=True
    )
    true_score, true_df_snp = python_score(
        plink_path, df_weight, weight_cols, center=True
    )
    assert np.allclose(df_score.values, true_score, atol=1e-4)
    assert df_score_snp.equals(true_df_snp)

    # apply to GBR data
    freq = pd.read_csv(plink_path.replace(".pgen", ".afreq"), sep="\t")[
        "ALT_FREQS"
    ].values

    # test --read-freq works properly when analyzing subset of the data.
    gbr_plink_path = dapgen.get_test_data("plink2.merged.GBR.pgen")
    shutil.copyfile(
        plink_path.replace(".pgen", ".afreq"), gbr_plink_path.replace(".pgen", ".afreq")
    )

    df_score, df_score_snp = dapgen.score(
        gbr_plink_path, df_weight, freq_suffix=".afreq", center=True
    )
    true_score, true_df_snp = python_score(
        gbr_plink_path, df_weight, weight_cols, freq=freq, center=True
    )
    os.remove(gbr_plink_path.replace(".pgen", ".afreq"))
    assert np.allclose(df_score.values, true_score, atol=1e-4)
    assert df_score_snp.equals(true_df_snp)


def test_two_snp():
    # build df_weight
    plink_path = dapgen.get_test_data("plink2.merged.pgen")
    np.random.seed(1234)
    df_weight = dapgen.read_pvar(plink_path.replace(".pgen", ".pvar"))
    weight_cols = ["WEIGHT1", "WEIGHT2", "WEIGHT3"]
    for col in weight_cols:
        df_weight[col] = np.random.normal(size=len(df_weight))

    # two individuals
    path1 = dapgen.get_test_data("plink2.chr21.pgen")
    path2 = dapgen.get_test_data("plink2.chr22.pgen")
    merged_path = dapgen.get_test_data("plink2.merged.pgen")
    plink_path1 = [path1, path2]
    plink_path2 = merged_path

    df_score, df_score_snp = dapgen.score(plink_path1, df_weight)
    true_score, true_df_snp = python_score(plink_path2, df_weight, weight_cols)
    assert np.allclose(df_score.values, true_score, rtol=1e-4, atol=1e-6)
    assert df_score_snp.equals(true_df_snp)
