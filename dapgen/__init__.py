__version__ = "0.1"

from ._read import read_pfile, read_pgen, read_pvar, read_psam, read_plink
from ._data import get_test_data
from ._utils import align_snp, score
from ._write import write_pgen

__all__ = [
    "__version__",
    "read_pfile",
    "read_pgen",
    "read_pvar",
    "read_psam",
    "align_snp",
    "get_test_data",
]
