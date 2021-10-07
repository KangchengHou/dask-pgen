From pandas-plink

```bash
for suffix in bed bim fam; do
    wget https://raw.githubusercontent.com/limix/pandas-plink/main/pandas_plink/test/data_files/data.${suffix} -O plink1.${suffix}
done
```

From plink-ng

```bash
wget https://github.com/chrchang/plink-ng/raw/master/2.0/Tests/TEST_PHASED_VCF/1kg_phase3_chr21_start.vcf.gz -O phased.vcf
plink2 --vcf phased.vcf --out phased
```