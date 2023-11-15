from distutils.version import LooseVersion
from io import open

import setuptools
from setuptools import setup

# Added support for environment markers in install_requires.
if LooseVersion(setuptools.__version__) < "36.2":
    raise ImportError("setuptools>=36.2 is required")

setup(
    name="dask-pgen",
    version="0.1",
    description="dask pgen",
    url="https://github.com/kangchenghou/dask-pgen",
    author="Kangcheng Hou",
    author_email="kangchenghou@gmail.com",
    packages=["dapgen"],
    setup_requires=["numpy", "Cython"],
    install_requires=[
        "dask",
        "pandas",
        "numpy",
        "natsort",
        "fire",
        "structlog",
        "Pgenlib"
        # "Pgenlib @ git+https://github.com/chrchang/plink-ng.git@fc08f4396c521ce48058804c627f8a3397b1e3f7#egg=pgenlib&subdirectory=2.0/Python",
    ],
    scripts=["bin/dapgen"],
)
