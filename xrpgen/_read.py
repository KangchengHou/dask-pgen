import pgenlib
import numpy as np
import xarray as xr


def read_pgen(path, snp_chunk=1024):

    from dask.array import concatenate, from_delayed
    from dask.delayed import delayed

    pgen = pgenlib.PgenReader(bytes(path, "utf8"))
    n_indiv = pgen.get_raw_sample_ct()
    n_snp = pgen.get_variant_ct()

    snp_chunk_xs = []
    snp_start = 0
    while snp_start < n_snp:
        snp_stop = min(snp_start + snp_chunk, n_snp)
        snp_chunk_xs.append(
            from_delayed(
                value=delayed(_read_pgen_chunk)(path, snp_start, snp_stop),
                shape=(n_indiv, n_snp, 2),
                dtype=np.float32,
            )
        )
        snp_start = snp_stop
    pgen.close()
    return concatenate(snp_chunk_xs, 1, True)


def _read_pgen_chunk(path, snp_start, snp_stop):
    with pgenlib.PgenReader(bytes(path, "utf8")) as pgen:
        n_indiv = pgen.get_raw_sample_ct()
        geno = np.empty([snp_stop - snp_start, n_indiv * 2], dtype=np.int32)
        pgen.read_alleles_range(snp_start, snp_stop, geno)

    return np.dstack([geno[:, ::2], geno[:, 1::2]]).swapaxes(0, 1).astype(np.float32)
