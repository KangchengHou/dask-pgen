# dask-pgen
Lazy loading PLINK genotype data (.pgen or .bed) file with dask-array.

# Installation
```bash
pip install -U cython numpy
git clone https://github.com/KangchengHou/dask-pgen.git
cd dask-pgen; pip install -e .
```

# Example
```bash
dapgen score \
    --plink <plink> \
    --weights <weights_path> \
    --out <out_path> \
    --center True # center the genotype or not, default is False
```