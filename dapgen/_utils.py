import pandas as pd
from typing import Tuple
import numpy as np


def align_snp(
    df1: pd.DataFrame, df2: pd.DataFrame
) -> Tuple[pd.Index, pd.Index, np.ndarray]:
    """
    Align two SNP dataframes by SNP, CHROM, POS, with the following 3 steps:
    1. match SNPs in `df1` and `df2` by `CHROM` and `POS`
    2. Try match REF and ALT columns in `df1` and `df2`, either
        REF1 = REF2, ALT1 = ALT2, or
        REF1 = ALT1, ALT2 = REF2
    3. {-1, 1} vectors indicating which location needs to be flipped will be also returned
    so that df1.loc[idx1, cols] and df2.loc[idx2, cols] * sign can be aligned.

    TODO: ambiguous alleles (A/T and C/G) will be removed.
    TODO: strand flip will be coped.

    Parameters
    ----------
    df1 : pd.DataFrame
        Dataframe 1 containing CHROM, POS, REF, ALT
    df2 : pd.DataFrame
        Dataframe 2 containing CHROM, POS, REF, ALT

    Returns
    -------
    df1_idx : pd.Index
        Index of df1
    df2_idx : pd.Index
        Index of df2
    flip_sign : np.ndarray
        Sign of the alignment
    """
    # check df1.index and df2.index are unique
    assert (
        df1.index.is_unique and df2.index.is_unique
    ), "df1.index and df2.index must be unique"
    required_cols = ["CHROM", "POS", "REF", "ALT"]
    # check required columns are in df1 and df2
    for df in [df1, df2]:
        assert set(required_cols).issubset(
            set(df.columns)
        ), f"df1 and df2 must both contain {','.join(required_cols)}"
    df1, df2 = df1[required_cols].copy(), df2[required_cols].copy()
    df1.index.name = "ID"
    df2.index.name = "ID"

    df_merged = pd.merge(
        df1.reset_index(),
        df2.reset_index(),
        on=["CHROM", "POS"],
        suffixes=["1", "2"],
    )
    noflip_index = (df_merged["REF1"] == df_merged["REF2"]) & (
        df_merged["ALT1"] == df_merged["ALT2"]
    )
    flip_index = (df_merged["REF1"] == df_merged["ALT2"]) & (
        df_merged["ALT1"] == df_merged["REF2"]
    )
    df_merged.loc[noflip_index, "flip_sign"] = 1
    df_merged.loc[flip_index, "flip_sign"] = -1

    df1_idx = pd.Index(df_merged.loc[noflip_index | flip_index, "ID1"].values)
    df2_idx = pd.Index(df_merged.loc[noflip_index | flip_index, "ID2"].values)
    flip_sign = df_merged.loc[noflip_index | flip_index, "flip_sign"].values

    return df1_idx, df2_idx, flip_sign
