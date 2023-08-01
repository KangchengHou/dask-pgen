import numpy as np
import pgenlib
import dask.array as da
import pandas as pd
from typing import List
from typing import Union


def write_pgen(path: str, geno: Union[da.Array, np.ndarray]):
    """
    Write a pgen file.

    Parameters
    ----------
    path : str
        Path to write the pgen file to.
    geno : da.Array
        Genotype data to write. geno either has shape
        (n_snp, n_indiv, 2) phasing information will be stored
        (n_snp, n_indiv) only genotype dosage will be stored
    """

    # replace na with -9
    geno[np.isnan(geno)] = -9
    if isinstance(geno, da.Array):
        geno = geno.compute()

    if geno.ndim == 3:
        n_snp, n_indiv = geno.shape[0:2]
        geno = geno.reshape(n_snp, n_indiv * 2)
        with pgenlib.PgenWriter(
            path.encode("utf-8"),
            n_indiv,
            n_snp,
            nonref_flags=False,
            hardcall_phase_present=True,
        ) as writer:
            writer.append_alleles_batch(geno.astype(np.int32), all_phased=True)
    elif geno.ndim == 2:
        n_snp, n_indiv = geno.shape[0:2]
        with pgenlib.PgenWriter(
            path.encode("utf-8"),
            n_indiv,
            n_snp,
            nonref_flags=False,
            hardcall_phase_present=False,
        ) as writer:
            writer.append_biallelic_batch(geno.astype(np.int8))


def write_pvar(path: str, df_pvar: pd.DataFrame, header: List[str] = None):
    """
    Write a pvar file. pvar must starts with #CHROM, snp index will be turned into ID,
    and have CHROM, POS, REF, ALT

    Parameters
    ----------
    path : str
        Path to write the pvar file to.
    pvar : pd.DataFrame
        pvar data to write. pvar has shape (n_snp, n_attr)
    header : List[str]
        Header of the pvar file. Each line must start with #
    """
    # write df_snp
    df_pvar = df_pvar.reset_index().rename(columns={"snp": "ID", "CHROM": "#CHROM"})

    FIXED_COLS = ["#CHROM", "POS", "ID", "REF", "ALT"]
    df_pvar = df_pvar[FIXED_COLS + [col for col in df_pvar.columns if col not in FIXED_COLS]]
    with open(path, "w") as f:
        if (header is not None) and len(header) > 0:
            assert np.all([h.startswith("#") for h in header])
            f.writelines("\n".join(header))
            f.write("\n")

    df_pvar.to_csv(path, sep="\t", index=False, mode="a", float_format="%.8g")