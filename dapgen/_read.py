import pgenlib
import numpy as np
from typing import Any, Dict
import pandas as pd

# TODO: instead of SNP chunk, use expected storage size of genotype


def read_pfile(pfile: str, phase=False, snp_chunk: int = 1024):
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

    """Read pgen

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
    geno : dask.array.core.Array (n_snp, n_indiv, 2) or (n_snp, n_indiv)
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


# TODO: support fast reading a transposed genotype matrix
# There may be currently a bug in pgenlib.pyx, which causes the transpose to
# fail.
def _read_pgen_chunk(
    path: str, snp_start: int, snp_stop: int, phase: bool, n_indiv: int = None
):
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
            geno = np.empty(
                [snp_stop - snp_start, pgen.get_raw_sample_ct() * 2], dtype=np.int32
            )
            pgen.read_alleles_range(snp_start, snp_stop, geno)
        else:
            geno = np.empty(
                [snp_stop - snp_start, pgen.get_raw_sample_ct()], dtype=np.int32
            )
            pgen.read_range(snp_start, snp_stop, geno)

    geno = np.ascontiguousarray(geno, np.float32)
    geno[geno == -9] = np.nan

    if phase:
        return geno.reshape(snp_stop - snp_start, n_indiv, 2)
    else:
        return geno


def read_pvar(path):
    skiprows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                break
            skiprows += 1

    df_pvar = (
        pd.read_csv(path, delim_whitespace=True, skiprows=skiprows)
        .rename(columns={"#CHROM": "CHROM", "ID": "snp"})
        .set_index("snp")
    )
    assert np.isin(
        ["CHROM", "POS", "REF", "ALT"], df_pvar.columns
    ).all(), "pvar file must have columns CHROM, POS, SNP, REF, and ALT"

    return df_pvar


def read_psam(path):
    skiprows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#IID"):
                break
            skiprows += 1

    df_psam = (
        pd.read_csv(path, delim_whitespace=True, skiprows=skiprows)
        .rename(columns={"#IID": "indiv"})
        .set_index("indiv")
    )

    return df_psam
