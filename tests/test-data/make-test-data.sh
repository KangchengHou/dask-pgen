# Make PLINK2 dataset.
# Toy data of first 99 individuals in the 1000 Genomes Project.
# plink2.chr21: randomly sampled SNPs in chr21
# plink2.chr22: randomly sampled SNPs in chr22
# plink2.merged: merged data of chr21 and chr22
# plink2.merged.GBR: GBR sub-sample of merged data (N=82)
# plink2.merged.FIN: FIN sub-sample of merged data (N=17)
# plink1.merged: merged data of chr21 and chr22, but in plink1 format

# chromosome 21
wget https://www.dropbox.com/s/ag53caz6s0kkcv1/chr21_phase3.pgen.zst?dl=1 -O raw.chr21.pgen.zst
wget https://www.dropbox.com/s/9xqkggh71guvup5/chr21_phase3_noannot.pvar.zst?dl=1 -O raw.chr21.pvar.zst
# chromosome 22
wget https://www.dropbox.com/s/ozraccaavbtdkzm/chr22_phase3.pgen.zst?dl=1 -O raw.chr22.pgen.zst
wget https://www.dropbox.com/s/g5sucurqv46y6q9/chr22_phase3_noannot.pvar.zst?dl=1 -O raw.chr22.pvar.zst
# sample file is shared for the two chromosomes
wget https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1 -O raw.psam

# extract the first 100 rows in column 1 of raw.psam
# (N = 99 because the first row is a header)
head -n 100 raw.psam | cut -f 1 >indiv.txt

for chrom in {21..22}; do
    plink2 --zst-decompress raw.chr"${chrom}".pgen.zst >raw.chr"${chrom}".pgen
    plink2 --zst-decompress raw.chr"${chrom}".pvar.zst >raw.chr"${chrom}".pvar
    cp raw.psam raw.chr"${chrom}".psam
    plink2 --pfile raw.chr"${chrom}" \
        --rm-dup exclude-all \
        --max-alleles 2 \
        --thin-count 2000 \
        --maf 0.01 \
        --keep indiv.txt \
        --snps-only \
        --seed 0 \
        --make-pgen --out plink2.chr"${chrom}"
done

# merge into one data set
plink2 --pfile plink2.chr21 --pmerge plink2.chr22 --out plink2.merged

# convert to plink1
plink2 --pfile plink2.merged --make-bed --out plink1.merged

# split by population
# extract the 1st column where the 4th column == GBR using awk, seperated by space
for pop in GBR FIN; do
    awk -v pop="$pop" '{if($4==pop) print $1}' plink2.merged.psam >indiv.txt
    plink2 --pfile plink2.merged \
        --keep indiv.txt \
        --make-pgen --out plink2.merged.${pop}
done

plink2 \
    --pfile plink2.merged \
    --freq \
    --out plink2.merged

# clean up
rm indiv.txt
rm *.zst
rm raw.*
rm *.log
