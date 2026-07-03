# Prepare input (compress + index)
bgzip populations.snps.vcf
bcftools index populations.snps.vcf.gz
# FIRST GET RID OF SAMPLES WITH CATASTROPHIC MISSINNESS < median_present_loci*0.1 (VERY LENIENT FILTER ON THE FULL ASSEMBLY)
stacks_dist_extract=/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/scripts/stacks-dist-extract
$stacks_dist_extract populations.log.distribs loci_per_sample > loci_per_sample.log.distribs
sort -nk1 loci_per_sample.log.distribs | awk '{
    sum += $3;
    vals[NR] = $3;
}
END {
    mean = sum / NR;
    if (NR % 2) {
        median = vals[(NR + 1) / 2];
    } else {
        median = (vals[NR / 2] + vals[NR / 2 + 1]) / 2;
    }
    printf "Mean: %.2f\nMedian: %.2f\n", mean, median;
}'
#cutoff:
#for n3 14686 * 0.3 = 4405
awk 'NR>1 && $3 < 4405 {print $1}' loci_per_sample.log.distribs > samples_to_remove.tmp.txt
# keep samples that are unique and interesting (despite having to much missing data)
cat samples_to_remove.tmp.txt | grep -v wasu_905_2 | grep -v judd_3639_1 | grep -v murp_930_1 > samples_removed.txt
#remove samples from assembly
bcftools view -S ^samples_removed.txt populations.snps.vcf.gz | bgzip -c > tmp1_raw_selSpec.snps.vcf.gz
bcftools index tmp1_raw_selSpec.snps.vcf.gz
#Filter 1:  
# Step 1: filter low depth genotypes
# Step 2: filter imbalanced het genotypes with binomial test
#     p-value of 0.01 tighten to 0.001 for more conservative (keep more borderline hets) 
#     or loosen to 0.05 if you want to be stricter about balance.
bcftools +setGT tmp1_raw_selSpec.snps.vcf.gz -- -t q -n . -i 'FMT/DP <= 3' | \
bcftools +setGT -- -t "b:AD<0.001" -n . | \
bgzip -c > filter1.tmp.vcf.gz
bcftools index filter1.tmp.vcf.gz
#Filter 2: SNPs with depth > 4× mean depth
#Compute mean depth:
bcftools query -f '[%DP\n]' filter1.tmp.vcf.gz | awk '$1 != "." {sum+=$1; count++} END {print sum/count}' > meanDP.tmp.txt
MEAN_DEPTH=$(cat meanDP.tmp.txt)
MAX_DEPTH=$(python3 -c "print($MEAN_DEPTH * 5)")
#for n3 mean depth is 20.75 (quite low since depth is variable i will use a much higher cutoff mean *10 = ~200)

MAX_DEPTH=200
#Filter sites:
bcftools +setGT filter1.tmp.vcf.gz \
  -- -t q -n . \
  -i "FMT/DP > $MAX_DEPTH" | bgzip -c  > filter2.tmp.vcf.gz
bcftools index filter2.tmp.vcf.gz
#ANOTATE HWE exact p-value in the INFO FIELD
bcftools +fill-tags filter2.tmp.vcf.gz -- -t F_MISSING,MAF,HWE | bgzip -c > filter2_HWEannotated.tmp.vcf.gz
bcftools index filter2_HWEannotated.tmp.vcf.gz
#Filter 3

#define parameters
MAX_F_MISSING=0.95
MAX_HWE=1e-15
MIN_MAF=0.005
prefix=SelSpec_FMISS095_MINMAC4

# Minor allele count >=4
# HWE > 1e-15 (very lenient) genus level, is just to filter spurious Alleles (eg, P allele freq ~0.5 with very low Het, not likely to be real)
# Remove SNPs missing in >95% of samples
#minimum samples per variant cutoff and prefix
bcftools view -i "MAF>=${MIN_MAF} && F_MISSING<=${MAX_F_MISSING} && HWE > ${MAX_HWE}" \
  filter2_HWEannotated.tmp.vcf.gz | bgzip -c > filter3.tmp.vcf.gz
bcftools index filter3.tmp.vcf.gz
#compute per sample and perloci stats
# COmpute and extract per site statistics into a table
bcftools +fill-tags filter3.tmp.vcf.gz -- -t F_MISSING,HWE,ExcHet,AC,AN,AF,MAF | 
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%F_MISSING\t%ExcHet\t%HWE\n' > site_stats.tmp.txt
#FILTER PARALOGS (LOCI WITH TOO MANY SNPS)
#count number of SNPs per loci
awk '{print $1}' site_stats.tmp.txt | sort | uniq -c | awk '{
    # Remove leading whitespace from uniq -c output
    gsub(/^[ \t]+/, "", $0)
    count = $1
    locus = $2
    print count, locus
}' > numberSNPs_per_loci.tmp.txt
sort -k1,1n numberSNPs_per_loci.tmp.txt | cut -d " " -f1 | uniq -c > ${prefix}_numberSNPs_per_loci_distribution.txt
#For rad seq > 4 (population level) is suspicious, but in this genus level assembly we many more.
#Distribution is not bimodal —> no clean trough separating "real loci" from "paralogs". 
#whitelist loci with less than 9 SNPs (the rest are likely paralog piled up) after which we see sharp decline
awk '$1 <= 8 {print $2}' numberSNPs_per_loci.tmp.txt > keep_loci.tmp.txt
#Since each contig spans positions 1–36, create a proper regions file:
awk '{print $1"\t1\t36"}' keep_loci.tmp.txt > regions.tmp.txt
#subset vcf by "contigs" (radloci)
bcftools view -R regions.tmp.txt filter3.tmp.vcf.gz | bgzip -c > ${prefix}.vcf.gz
bcftools index ${prefix}.vcf.gz
#check how many variants remain
bcftools query -f '%CHROM' ${prefix}.vcf.gz | wc -l
#check how many loci remain
bcftools query -f '%CHROM' ${prefix}.vcf.gz | sort | uniq | wc -l


#Now also make a dataset keeping only one snp per loci
#keep one snp per loci "LD prune" bcftools +prune -w 36 -m 1  is not working
#so i will keep one snp per loci (the one with the highest AF)
# Extract CHROM, POS, and AF, then sort by AF (highest first)
bcftools query -f '%CHROM\t%POS\t%AF\n'  ${prefix}.vcf.gz | sort -k3,3 -nr > ${prefix}.tmp1.AFsorted_snps.txt
# 4. Keep only the first occurrence of each chromosome = the one with highest AF
awk '!seen[$1]++' ${prefix}.tmp1.AFsorted_snps.txt | awk '{print $1"\t"$2"\t"$2}' > ${prefix}.tmp1.keep_snps.txt

#keep only those snps
bcftools view -R ${prefix}.tmp1.keep_snps.txt ${prefix}.vcf.gz | bgzip -c > ${prefix}.usnps.vcf.gz
bcftools index ${prefix}.usnps.vcf.gz

#count how many snps left, should also match number of unique CHROM
bcftools index -n ${prefix}.usnps.vcf.gz
bcftools query -f '%CHROM\n' ${prefix}.usnps.vcf.gz | uniq | wc -l

#CONVERT TO PHYLIP USING VCF2PHYLIP
#bcftools query -l ${prefix}_selSpec.gz | awk '{print $1 "\t" $1}' > ${prefix}_selSpec_individual_popmap.txt
#vcf2phylip=~/tools/vcf2phylip/vcf2phylip.py
#python3 $vcf2phylip -i ${prefix}_selSpec.vcf.gz

vcf2phylip=/storage1/fs1/christine.e.edwards/Active/aguirre/vcf2phylip/vcf2phylip.py
python3 $vcf2phylip -i min50_selSpec.vcf.gz
python3 $vcf2phylip -i min40_selSpec.vcf.gz
python3 $vcf2phylip -i min30_selSpec.vcf.gz

#download from cluster all assembly outputs
scp -r a.eduardo@compute1-client-1.ris.wustl.edu:/storage1/fs1/christine.e.edwards/Active/weinmannia/assembly/populations .
