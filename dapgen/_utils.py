import pandas as pd
from typing import Tuple, List, Optional
import numpy as np
import tempfile
import os
import subprocess
import shutil
import os
from contextlib import contextmanager
import os
from tqdm import tqdm
from ._read import read_pvar, read_bim, parse_plink_path, infer_merge_dim


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def get_dependency(name, download=True):
    """Get path to an depenency
    Find the binary in the following locations:
    - $PATH
    - package installment directory admix-tools/.admix_cache/bin/<name>
    If not found in any of these locations, download the corresponding software package
    - plink: https://www.cog-genomics.org/plink/2.0/
    Parameters
    ----------
    download : bool
        whether to download plink if not found
    Returns
    -------
    Path to binary executable
    """

    # find in path
    if shutil.which(name):
        return shutil.which(name)
    else:
        raise ValueError(
            f"{name} not found in PATH. "
            "Download plink2 (https://www.cog-genomics.org/plink/2.0/) "
            "and put it in the PATH"
        )


def _score_single_plink(
    path: str,
    df_weight: pd.DataFrame,
    weight_cols: List[str],
    center: bool,
    freq_suffix: Optional[str],
    **kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Score a single plink file against a pd.DataFrame

    Parameters
    ----------
    path : str
        [description]
    df_weight : pd.DataFrame
        [description]

    Notes
    -----

    In some cases (e.g., All of Us), we have two splitted SNPs in genotype file
    chr1:51001424:C:CT and chr1:51001424:CT:C, this function will match both SNPs using
    the same row from weight file. `df_weight` will contain two rows
    (chr1:51001424:C:CT, chr1:51001424:CT:C) with opposite effects, just as if there
    are indeed two SNPs. Then everything else would follow through.
    """
    if path.endswith(".pgen"):
        df_snp = read_pvar(path.replace(".pgen", ".pvar"))
    elif path.endswith(".bed"):
        df_snp = read_bim(path.replace(".bed", ".bim"))

    # using CHROM, POS, REF, ALT to align
    # flip_sign is not used because plink2 --score will handle this
    df_snp_idx, df_weight_idx, flip_sign = align_snp(df_snp, df_weight)

    # subset matched snps, since df_weight will be subsetted, plink2 will cope filter
    # the `df_snp`, so we don't need to subset `df_snp` here
    df_snp = df_snp.loc[df_snp_idx]
    df_weight = df_weight.loc[df_weight_idx, ["ALT"] + weight_cols]
    df_weight.index = df_snp_idx
    df_weight.index.name = "SNP"

    plink2_bin = get_dependency("plink2")
    with tempfile.TemporaryDirectory() as tmp_dir:
        df_weight.to_csv(os.path.join(tmp_dir, "weight.txt"), sep="\t", index=True)
        np.savetxt(
            os.path.join(tmp_dir, "weight_snps.txt"),
            df_snp_idx,
            fmt="%s",
            delimiter="\n",
        )
        cmds = [
            f"{plink2_bin} --extract {os.path.join(tmp_dir, 'weight_snps.txt')} --score {tmp_dir}/weight.txt 1 2 header-read cols=+scoresums,-scoreavgs"
        ]
        if center:
            # append center to --score argument
            cmds[-1] += " center"

        # keyword arguments
        add_cmds = [f" --{k.replace('_', '-')} {kwargs[k]}" for k in kwargs]
        cmds += add_cmds

        max_weight_col = len(df_weight.columns) + 1
        if max_weight_col > 3:
            cmds += [f"--score-col-nums 3-{max_weight_col}"]
        else:
            cmds += [f"--score-col-nums 3"]

        if path.endswith(".pgen"):
            prefix = path[:-5]
            cmds += [f"--pfile {prefix}"]
        elif path.endswith(".bed"):
            prefix = path[:-4]
            cmds += [f"--bfile {path[:-4]}"]

        # specify frequency file if available
        if freq_suffix is not None:
            cmds += [f"--read-freq {prefix + freq_suffix}"]

        cmds += [f"--out {tmp_dir}/out"]

        subprocess.check_call(" ".join(cmds), shell=True)

        # read back in
        df_score = pd.read_csv(os.path.join(tmp_dir, "out.sscore"), sep="\t")
    if path.endswith(".pgen"):
        df_score = df_score.set_index(df_score.columns[0])
        df_score.index.name = "indiv"
    elif path.endswith(".bed"):
        df_score.index = (
            df_score.iloc[:, 0].astype(str) + "_" + df_score.iloc[:, 1].astype(str)
        )
        df_score.index.name = "indiv"
    else:
        raise ValueError("path must end with .pgen or .bed")
    df_score = df_score.loc[:, [col + "_SUM" for col in weight_cols]]
    df_score.columns = weight_cols
    return df_score, df_snp


def freq(plink_path: str, memory: Optional[int] = None) -> pd.DataFrame:
    """Calculate the frequency

    Parameters
    ----------
    plink_path : str
        path to the plink file
    memory: int
        memory constraint.

    Returns
    -------
    pd.DataFrame
        with columns #CHROM ID REF ALT ALT_FREQS OBS_CT
    """
    plink2_bin = get_dependency("plink2")
    with tempfile.TemporaryDirectory() as tmp_dir:
        cmds = [f"{plink2_bin} --freq --out {tmp_dir}/freq"]
        if memory is not None:
            cmds += [f"--memory {memory * 1024}"]

        if plink_path.endswith(".pgen"):
            prefix = plink_path[:-5]
            cmds += [f"--pfile {prefix}"]
        elif plink_path.endswith(".bed"):
            prefix = plink_path[:-4]
            cmds += [f"--bfile {plink_path[:-4]}"]

        subprocess.check_call(" ".join(cmds), shell=True)

        # read back in
        df_freq = pd.read_csv(os.path.join(tmp_dir, "freq.afreq"), sep="\t")

    return df_freq


def _score_multiple_plink(
    path_list: List[str],
    df_weight: pd.DataFrame,
    weight_cols: List[str],
    center: bool,
    freq_suffix: Optional[str],
    **kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # multiple plink files
    merge_dim = infer_merge_dim(path_list)
    if merge_dim == "indiv":
        raise ValueError("indiv merge not supported")

    assert merge_dim in ["snp", "indiv"]

    if merge_dim == "indiv":
        df_score_list = []
        df_snp = None
        for path in tqdm(path_list, desc="scoring multiple plink files"):
            df_score, this_df_snp = _score_single_plink(
                path=path,
                df_weight=df_weight,
                weight_cols=weight_cols,
                center=center,
                freq_suffix=freq_suffix,
                **kwargs,
            )
            df_score_list.append(df_score)
            if df_snp is None:
                df_snp = this_df_snp
            else:
                assert df_snp.equals(
                    this_df_snp
                ), "df_snp must be the same for all indiv"
    elif merge_dim == "snp":
        df_score = None
        df_snp_list = []
        for path in tqdm(path_list, desc="scoring multiple plink files"):
            this_df_score, df_snp = _score_single_plink(
                path=path,
                df_weight=df_weight,
                weight_cols=weight_cols,
                center=center,
                freq_suffix=freq_suffix,
                **kwargs,
            )
            if df_score is None:
                df_score = this_df_score
            else:
                assert df_score.index.equals(
                    this_df_score.index
                ), "index must be the same"
                df_score += this_df_score
            df_snp_list.append(df_snp)
        df_snp = pd.concat(df_snp_list, axis=0)
    else:
        raise ValueError("merge_dim must be 'indiv' or 'snp'")
    return df_score, df_snp


def score(
    plink_path: str,
    df_weights: pd.DataFrame,
    weight_cols: List[str] = None,
    center: bool = False,
    freq_suffix: Optional[str] = None,
    **kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Wrapper for scoring plink files against a pd.DataFrame

    CHROM, POS, REF, ALT columns must be present in df_weights and will be used to align
    with the plink file.

    All other numerical columns will be used as input for scoring.

    Parameters
    ----------
    plink_path : str
        Path to plink file
    df_weight : pd.DataFrame
        Dataframe containing weights
    out : str
        Path to output file
    freq_suffix : str
        if provided, every plink file [plink] must be associated with a [plink].[freq_suffix] file
    """
    # basic checks on df_weight
    assert np.all([col in df_weights.columns for col in ["CHROM", "POS", "REF", "ALT"]])

    if weight_cols is None:
        # get numerical columns except for CHROM, POS, REF, ALT
        weight_cols = [
            col
            for col in df_weights.select_dtypes(include=np.number).columns
            if col not in ["CHROM", "POS", "REF", "ALT"]
        ]

    # parse plink path
    path_list = parse_plink_path(plink_path)
    if len(path_list) == 1:
        # single plink file
        df_score, df_snp = _score_single_plink(
            path_list[0],
            df_weights,
            weight_cols,
            center=center,
            freq_suffix=freq_suffix,
            **kwargs,
        )
    elif len(path_list) > 1:
        # multiple plink files
        df_score, df_snp = _score_multiple_plink(
            path_list,
            df_weights,
            weight_cols,
            center=center,
            freq_suffix=freq_suffix,
            **kwargs,
        )
    return df_score, df_snp


def align_snp(
    df1: pd.DataFrame, df2: pd.DataFrame, by: str = "pos"
) -> Tuple[pd.Index, pd.Index, np.ndarray]:
    """
    Align two SNP dataframes by SNP, CHROM, POS, with the following 3 steps:
    1. match SNPs in `df1` and `df2` by `CHROM` and `POS`
    2. Try match REF and ALT columns in `df1` and `df2`, either
        REF1 = REF2, ALT1 = ALT2, or
        REF1 = ALT1, ALT2 = REF2
    3. {-1, 1} vectors indicating which location needs to be flipped will be also returned
    so that df1.loc[idx1, cols] and df2.loc[idx2, cols] * sign can be aligned.

    Parameters
    ----------
    df1 : pd.DataFrame
        Dataframe 1 containing CHROM, POS, REF, ALT
    df2 : pd.DataFrame
        Dataframe 2 containing CHROM, POS, REF, ALT
    by : str
        either by pos (chromosome and position) of snp (SNP ID)
    Returns
    -------
    df1_idx : pd.Index
        Index of df1
    df2_idx : pd.Index
        Index of df2
    flip_sign : np.ndarray
        Sign of the alignment
    """
    if not df1.index.is_unique:
        # print the first two row numbers that are not unique
        dup_index = df1.index[df1.index.duplicated(keep=False)]
        raise ValueError(
            f"df1.index must be unique, but contains {len(dup_index)} duplicates (e.g., {dup_index[0]})"
        )
    if not df2.index.is_unique:
        dup_index = df2.index[df2.index.duplicated(keep=False)]
        raise ValueError(
            f"df2.index must be unique, but contains {len(dup_index)} duplicates (e.g., {dup_index[0]})"
        )

    if by == "pos":
        required_cols = ["CHROM", "POS", "REF", "ALT"]
        on_cols = ["CHROM", "POS"]
    elif by == "id":
        required_cols = ["SNP", "REF", "ALT"]
        on_cols = ["SNP"]
    else:
        raise ValueError("by must be either 'pos' or 'id'")
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
        on=on_cols,
        suffixes=["1", "2"],
    )
    assert len(df_merged) > 0, "df1 and df2 must contain at least one common SNP"
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
