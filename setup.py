from distutils.version import LooseVersion
from io import open

import setuptools
from setuptools import setup

# Added support for environment markers in install_requires.
if LooseVersion(setuptools.__version__) < "36.2":
    raise ImportError("setuptools>=36.2 is required")

setup(
    name="xarray-pgen",
    version="0.1",
    description="xarray pgen",
    url="https://github.com/kangchenghou/xarray-pgen",
    author="Kangcheng Hou",
    author_email="kangchenghou@gmail.com",
    packages=["xrpgen"],
    setup_requires=["numpy"],
)
