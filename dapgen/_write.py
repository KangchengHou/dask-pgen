import numpy as np
import pgenlib
import dask.array as da


def write_pgen(path: str, geno: da.Array):
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
            writer.append_alleles_batch(
                geno.astype(np.int32).compute(), all_phased=True
            )
    elif geno.ndim == 2:
        n_snp, n_indiv = geno.shape[0:2]
        with pgenlib.PgenWriter(
            path.encode("utf-8"),
            n_indiv,
            n_snp,
            nonref_flags=False,
            hardcall_phase_present=False,
        ) as writer:
            writer.append_biallelic_batch(geno.astype(np.int8).compute())
