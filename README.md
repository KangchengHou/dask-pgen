# dask-pgen
Lazy loading PLINK genotype data (.pgen or .bed) file with dask-array.

# Installation
```bash
pip install -U cython numpy
git clone https://github.com/KangchengHou/dask-pgen.git
cd dask-pgen; pip install -e .
```
Or,
```bash
pip install git+https://github.com/KangchengHou/dask-pgen.git
```

# Example
```bash
# if dapgen is not found, set the path of dapgen yourself
# chmod +x dask-pgen/bin/dapgen
# dapgen=dask-pgen/bin/dapgen
# replace "dapgen" with "$dapgen"

dapgen score \
    --plink <plink> \
    --weights <weights_path> \
    --weight-col-prefix <weight_column_prefix> \
    --chrom-col CHR --pos-col POS --alt-col A1 --ref-col A2 \
    --out <out_path> \
    --center True # center the genotype or not, default is False \
    --threads 4 \
    --memory 20000
```
