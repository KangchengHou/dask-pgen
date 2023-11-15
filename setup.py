from distutils.version import LooseVersion
from io import open

from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="dask-pgen",
    version="0.1.1",
    description="Dask interface for memory-efficient reading PLINK genotype files",
    long_description=long_description,
    long_description_content_type="text/markdown",
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
        "Pgenlib",
    ],
    scripts=["bin/dapgen"],
)
