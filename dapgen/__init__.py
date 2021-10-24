__version__ = "0.1"

from ._read import read_pfile, read_pgen, read_pvar, read_psam
from ._data import get_test_data

__all__ = [
    "__version__",
    "read_pfile",
    "read_pgen",
    "read_pvar",
    "read_psam",
    "get_test_data",
]
