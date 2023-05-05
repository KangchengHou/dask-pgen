import pgenlib
import numpy as np
from typing import Any, Dict, Tuple, List, Union, Optional
import pandas as pd
import dask.array as da
import glob
import os
from natsort import natsorted


def infer_merge_dim(paths: List[str]):
    """
    Infer the dimension to merge for multiple plink files.
    """
    # infer merge_dim from the first two files
    if paths[0].endswith(".pgen"):
        df_snp1 = read_pvar(paths[0].replace(".pgen", ".pvar"))
        df_snp2 = read_pvar(paths[1].replace(".pgen", ".pvar"))
    elif paths[0].endswith(".bed"):
        df_snp1 = read_bim(paths[0].replace(".bed", ".bim"))
        df_snp2 = read_bim(paths[1].replace(".bed", ".bim"))
    else:
        raise ValueError("path must end with .pgen or .bed")
    if df_snp1.index.equals(df_snp2.index):
        merge_dim = "indiv"
    else:
        merge_dim = "snp"
    return merge_dim


def _read_single_plink(path: str, phase: bool, snp_chunk: int):
    geno = read_pgen(
        path,
        phase=phase,
        snp_chunk=snp_chunk,
    )
    if path.endswith(".pgen"):
        df_snp = read_pvar(path.replace(".pgen", ".pvar"))
        df_indiv = read_psam(path.replace(".pgen", ".psam"))
    elif path.endswith(".bed"):
        df_snp = read_bim(path.replace(".bed", ".bim"))
        df_indiv = read_fam(path.replace(".bed", ".fam"))
    return geno, df_snp, df_indiv


def _read_multiple_plink(
    paths: List[str], phase: bool, snp_chunk: int, merge_dim: Optional[str] = None
) -> Tuple[da.Array, pd.DataFrame, pd.DataFrame]:
    """_read_multiple_plink

    Parameters
    ----------
    paths : List[str]
        list of paths to plink files
    phase : bool
        Whether to read the phasing information of data
    snp_chunk : int
        number of SNPs within a chunk
    merge_dim : bool, optional
        dimension over which to merge, by default None, can be "snp" or "indiv"
        merge_dim = None: infer from the first two files
        merge_dim = "snp" will concatenate different SNP set into one array
        merge_dim = "indiv" will concatenate different individual set into one array

    Returns
    -------
    Tuple[da.Array, pd.DataFrame, pd.DataFrame]
        tuple of genotype, df_snp and df_indiv

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    """
    assert len(paths) > 1, "Need at least two files to merge"
    if merge_dim is None:
        merge_dim = infer_merge_dim(paths)

    assert merge_dim in ["snp", "indiv"]

    geno_list = []
    df_snp_list = []
    df_indiv_list = []
    for path in paths:
        geno, df_snp, df_indiv = _read_single_plink(path=path, phase=phase, snp_chunk=snp_chunk)
        geno_list.append(geno)
        df_snp_list.append(df_snp)
        df_indiv_list.append(df_indiv)

    if merge_dim == "indiv":
        # check snp are all the same
        assert np.all(
            [df_snp.equals(df_snp_list[0]) for df_snp in df_snp_list[1:]]
        ), "df_snp are not all the same when merging on indiv"
        df_snp = df_snp_list[0]
        geno = da.concatenate(geno_list, axis=1)
        df_indiv = pd.concat(df_indiv_list, axis=0)
    elif merge_dim == "snp":
        # determine the df_snp order such that the order is sorted by
        # last element of CHROM and POS in each df_snp
        df_chrom_pos = pd.DataFrame(
            [[df_snp.iloc[-1]["CHROM"], df_snp.iloc[-1]["POS"]] for df_snp in df_snp_list],
            columns=["CHROM", "POS"],
        )
        # sort by CHROM and POS
        order = df_chrom_pos.sort_values(by=["CHROM", "POS"]).index
        # adjust order for geno_list, df_snp_list
        geno_list = [geno_list[i] for i in order]
        df_snp_list = [df_snp_list[i] for i in order]
        assert np.all(
            [df_indiv.equals(df_indiv_list[0]) for df_indiv in df_indiv_list[1:]]
        ), "df_indiv are not all the same when merging on snp"
        # concatenate
        geno = da.concatenate(geno_list, axis=0)
        df_snp = pd.concat(df_snp_list, axis=0)
        df_indiv = df_indiv_list[0]
    else:
        raise ValueError("merge_dim must be either 'snp' or 'indiv'")

    # check df_snp is sorted by CHROM and POS
    assert df_snp.sort_values(by=["CHROM", "POS"]).equals(
        df_snp
    ), "df_snp is not sorted by CHROM and POS"
    return geno, df_snp, df_indiv


def parse_plink_path(pathname: Union[str, List]) -> List[str]:
    """Parse a input plink file pattern into a list of files

    Parameters
    ----------
    filename : str
        Path to the plink file

    Returns
    -------
    list
        List of lists containing the plink file
    """
    if isinstance(pathname, list):
        # already a list of paths
        out = pathname
    elif isinstance(pathname, str):
        if ("*" in pathname) or os.path.isdir(pathname):
            # find all files matching the pattern
            if "*" in pathname:
                # pattern is a file pattern
                out = glob.glob(pathname)
            elif os.path.isdir(pathname):
                # pathname is a directory
                out = out = glob.glob(pathname + "/*")

            pgen_list = [p for p in out if p.endswith(".pgen")]
            bed_list = [p for p in out if p.endswith(".bed")]
            assert (len(pgen_list) > 0) != (
                len(bed_list) > 0
            ), f"Either the directory={pathname} contains .pgen or .bed files, but not both"
            if len(pgen_list) > 0:
                out = pgen_list
            elif len(bed_list) > 0:
                out = bed_list
            else:
                raise ValueError(
                    f"Either the directory={pathname} contains .pgen or .bed files," " but not both"
                )
            out = natsorted(out)
        elif pathname.endswith(".bed") or pathname.endswith(".pgen"):
            out = [pathname]
        elif pathname.endswith(".txt"):
            # pattern is a file with a list of files
            with open(pathname, "r") as f:
                out = [line.strip() for line in f]
        else:
            raise ValueError("Unable to parse plink pathname")

    # either all path endswith .pgen or .bed
    assert all(p.endswith(".pgen") for p in out) or all(
        p.endswith(".bed") for p in out
    ), "Either all files end with .pgen or all files end with .bed"
    return out


def read_plink(pathname, phase=False, snp_chunk=1024):
    """General-purpose function to read plink files
    Usage includes
    - read_plink(pathname="chr21.bed", phase=False)
    - read_plink(pathname="chr22.pgen", phase=True)
    - read_plink(pathname="*.pgen", phase=True)
    - read_plink(pathname="file_list.txt") # file_list.txt contains rows of file names

    Parameters
    ----------
    pathname : str
        Path to plink file prefix without .pgen/.pvar/.psam
    phase : bool
        Whether to read the phasing information of data
    snp_chunk: int
        number of SNPs within a chunk

    Returns
    -------
    (geno, df_snp, df_indiv)
        geno: genotype matrix
            (n_snp, n_indiv) if phase is set to False
            (n_snp, n_indiv, 2) if phase is set to True
        df_snp: SNP information data frame
        df_indiv: individual information data frame
    """

    path_list = parse_plink_path(pathname)
    if len(path_list) == 1:
        # only one file
        return _read_single_plink(path_list[0], phase=phase, snp_chunk=snp_chunk)
    else:
        # multiple files
        return _read_multiple_plink(path_list, phase=phase, snp_chunk=snp_chunk)


# TODO: after read_plink is done, remove this function
# TODO: instead of SNP chunk, use expected storage size of genotype
def read_pfile(
    pfile: str, phase=False, snp_chunk: int = 1024
) -> Tuple[da.Array, pd.DataFrame, pd.DataFrame]:
    """read plink file and form xarray.Dataset

    Parameters
    ----------
    pfile : str
        path to plink file prefix without .pgen/.pvar/.psam
    phase : bool
        whether to read the phasing information of data
    snp_chunk: int
        number of SNPs within a chunk

    Returns
    -------
    (geno, df_snp, df_indiv)
        geno: genotype matrix
            (n_snp, n_indiv) if phase is set to False
            (n_snp, n_indiv, 2) if phase is set to True
        df_snp: SNP information data frame
        df_indiv: individual information data frame
    """

    # count number of a0 as dosage, (A1 in usual PLINK bim file)
    geno = read_pgen(
        pfile + ".pgen",
        phase=phase,
        snp_chunk=snp_chunk,
    )
    df_snp = read_pvar(pfile + ".pvar")
    df_indiv = read_psam(pfile + ".psam")

    return geno, df_snp, df_indiv


def read_pgen(path: str, snp_chunk: int = 1024, phase: bool = False):

    """Read pgen file including .pgen and .bed files

    Parameters
    ----------
    path : str
        Path to pgen file, e.g. /path/to/file.pgen or /path/to/file.bed
        suffix needs to be .pgen or .bed
    snp_chunk : int
        Number of SNPs to read at once
    phase : bool
        If True, return phased genotypes (n_indiv, n_snp, 2)
        If False, return unphased genotypes (n_indiv, n_snp)

    Returns
    -------
    geno : da.Array (n_snp, n_indiv, 2) or (n_snp, n_indiv)
    """

    from dask.array import concatenate, from_delayed
    from dask.delayed import delayed

    # raw_sample_ct is required for a .bed file, and optional for a .pgen file
    if path.endswith(".bed"):
        n_indiv: int = sum(1 for line in open(path.replace(".bed", ".fam")))
        pgen = pgenlib.PgenReader(bytes(path, "utf8"), raw_sample_ct=n_indiv)
    else:
        pgen = pgenlib.PgenReader(bytes(path, "utf8"), raw_sample_ct=None)

    n_indiv = pgen.get_raw_sample_ct()
    n_snp = pgen.get_variant_ct()

    snp_chunk_xs = []
    snp_start = 0
    while snp_start < n_snp:
        snp_stop = min(snp_start + snp_chunk, n_snp)
        if phase:
            shape: Any = (snp_stop - snp_start, n_indiv, 2)
        else:
            shape = (snp_stop - snp_start, n_indiv)

        snp_chunk_xs.append(
            from_delayed(
                value=delayed(_read_pgen_chunk)(
                    path=path,
                    snp_start=snp_start,
                    snp_stop=snp_stop,
                    phase=phase,
                    n_indiv=n_indiv,
                ),
                shape=shape,
                dtype=np.float32,
            )
        )
        snp_start = snp_stop
    pgen.close()
    return concatenate(snp_chunk_xs, 0, False)


def _read_pgen_chunk(path: str, snp_start: int, snp_stop: int, phase: bool, n_indiv: int = None):
    """
    Read a chunk of SNPs from a pgen file

    Parameters
    ----------
    path : str
        Path to pgen file, e.g. /path/to/file.pgen or /path/to/file.bed
        suffix needs to be .pgen or .bed
    snp_start : int
        First SNP to read
    snp_stop : int
        Last SNP to read
    phase : bool
        If True, return phased genotypes (n_snp, n_indiv, 2)
        If False, return unphased genotypes (n_snp, n_indiv)
    n_indiv : int
        Number of individuals in the pgen file
    """
    # np.int32 is required by pgenlib
    with pgenlib.PgenReader(bytes(path, "utf8"), raw_sample_ct=n_indiv) as pgen:
        if phase:
            geno = np.empty([snp_stop - snp_start, pgen.get_raw_sample_ct() * 2], dtype=np.int32)
            pgen.read_alleles_range(snp_start, snp_stop, geno)
        else:
            geno = np.empty([snp_stop - snp_start, pgen.get_raw_sample_ct()], dtype=np.int32)
            pgen.read_range(snp_start, snp_stop, geno)

    geno = np.ascontiguousarray(geno, np.float32)
    geno[geno == -9] = np.nan

    if phase:
        return geno.reshape(snp_stop - snp_start, n_indiv, 2)
    else:
        return geno


def read_bim(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        header=None,
        delim_whitespace=True,
        names=["CHROM", "snp", "CM", "POS", "ALT", "REF"],
    ).set_index("snp")


def read_fam(path: str):
    df = pd.read_csv(
        path,
        header=None,
        delim_whitespace=True,
        usecols=[0, 1],
        names=["FID", "IID"],
    ).astype(str)
    df.index = df.iloc[:, 0].astype(str) + "_" + df.iloc[:, 1].astype(str)
    df.index.name = "indiv"
    return df


def read_pvar(path: str, return_header: bool = False):
    """
    Read pvar file

    Parameters
    ----------
    path : str
        Path to pvar file, e.g. /path/to/file.pvar
    return_header : bool
        If True, return header (commented line starting with #) as a list of strings
    """
    header = []
    skiprows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                break
            header.append(line.strip())
            skiprows += 1

    df_pvar = (
        pd.read_csv(path, delim_whitespace=True, skiprows=skiprows)
        .rename(columns={"#CHROM": "CHROM", "ID": "snp"})
        .set_index("snp")
    )
    assert np.isin(
        ["CHROM", "POS", "REF", "ALT"], df_pvar.columns
    ).all(), "pvar file must have columns CHROM, POS, SNP, REF, and ALT"

    if return_header:
        return df_pvar, header
    else:
        return df_pvar


def read_psam(path):
    skiprows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#IID") or line.startswith("#FID"):
                break
            skiprows += 1

    df_psam = pd.read_csv(path, delim_whitespace=True, skiprows=skiprows)
    # only one column as index
    if df_psam.columns[0] == "#IID":
        df_psam = df_psam.set_index("#IID")
    elif df_psam.columns[0] == "#FID":
        assert (
            df_psam.columns[1] == "IID"
        ), "if psam's first column is #FID, the second column must be IID"
        df_psam.index = df_psam["#FID"].astype(str) + "_" + df_psam["IID"].astype(str)
        df_psam = df_psam.drop(columns=["#FID", "IID"])
    df_psam.index = df_psam.index.astype(str)
    df_psam.index.name = "indiv"
    return df_psam
