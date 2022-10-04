__version__ = "0.1"

from ._read import (
    read_pfile,
    read_pgen,
    read_pvar,
    read_psam,
    read_bim,
    read_fam,
    read_plink,
    parse_plink_path,
)
from ._data import get_test_data
from ._utils import align_snp, score, freq
from ._write import write_pgen, write_pvar

__all__ = [
    "__version__",
    "read_pfile",
    "read_pgen",
    "read_pvar",
    "read_psam",
    "write_pgen",
    "write_pvar",
    "align_snp",
    "get_test_data",
]
