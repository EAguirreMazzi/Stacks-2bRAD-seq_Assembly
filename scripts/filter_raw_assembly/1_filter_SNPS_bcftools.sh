#NEW FILTERING APPROACH:
#this pipeline is to match out divide and conquer approach
#the isdea is to filter low quality-spurious loci from the full assembly 
#no metric accounting for allele balance included
#this filtered dataset will be used to prepare subset datasets for regions 
#(individual level filters and AF filters will be applied within Subsets)
#define parameters
#!/bin/bash
input_vcf=$1
prefix="$2"
#1st remove samples with catastrophic missing data (since they will bias cutoff downstream)
# REMOVE SAMPLES WITH CATASTROPHIC MISSINNESS > 98% (VERY LENIENT FILTER)
# compute PSC stats table
nSNPs1=$(bcftools query -f '%CHROM' ${input_vcf} | wc -l)
#compute missingness per sample
bcftools stats -s - ${input_vcf} \
| grep "^PSC" \
| awk -v n=$nSNPs1 '{pMissing=$14 / n; print $3 "\t" pMissing}' \
> temp0.${prefix}.pMissing
#explore missingnes distribution (EXPLORE TO LEARN ABOUT THE DATA and decide cutoff)
awk '{dist[sprintf("%.2f", $2)]++} END {for (n in dist) print dist[n], n}' temp0.${prefix}.pMissing \
| sort -k2,2n \
> ${prefix}_raw.pMissing.distr
#FLAG specimens with to much missing data and remove from dataset
awk '$2 > 0.99 {print $1}' temp0.${prefix}.pMissing > ${prefix}.temp0.REMOVED_SAMPLES.txt
# keep samples that are unique and interesting (despite having to much missing data, for phylogeny only)
cat ${prefix}.temp0.REMOVED_SAMPLES.txt | grep -v wasu_905_2 | grep -v judd_3639_1 | grep -v murp_930_1 > ${prefix}_REMOVED_SAMPLES.txt
#apply the filter
bcftools view -s ^$(paste -sd, ${prefix}_REMOVED_SAMPLES.txt) ${input_vcf} | bcftools view -i "MAC>1" | bgzip -c >  ${prefix}.temp0.vcf.gz
bcftools index ${prefix}.temp0.vcf.gz
#NOW remove sites with low depth < 3 (unreliable genotype calling) and too high depth (likely paralogs stacked)
#ANNOTATE (or UPDATE) TOTAL DEPTH PER SNP and NS (number of samples)
bcftools +fill-tags ${prefix}.temp0.vcf.gz -- -t 'DP:1=int(sum(FORMAT/DP)),NS' | bgzip -c > ${prefix}.temp1.vcf.gz
bcftools index ${prefix}.temp1.vcf.gz
#COMPUTE NORMALIZED DEPTH (DP/NS)
#Explore Numerical Distribution of DP to decide cutoff
bcftools query -f '%INFO/DP\t%INFO/NS\n' ${prefix}.temp1.vcf.gz \
| awk '{a[int($1/$2)]++} END {for (i in a) print a[i], i}' \
| sort -k2,2n \
> ${prefix}_SNP_depth_distribution.txt
#Compute the maxdepth cutof f= weighted mean percentile 99th:
MAX_DEPTH=$(awk '
{
  sum += $1*$2; n += $1
  for (i=1; i<=$1; i++) print $2
}' ${prefix}_SNP_depth_distribution.txt \
| sort -n \
| awk '
BEGIN { total=0 }
{ vals[++total] = $1 }
END { print vals[int(total*0.99)] }')
echo "MAX_DEPTH (99th percentile): $MAX_DEPTH"
#alternatively use mean depth*3
#MAX_DEPTH=$(awk '{sum += $1*$2; n += $1} END {print (sum/n)*3}' ${prefix}_SNP_depth_distribution.txt)
#echo "maxdepth cutoff= weighted mean depth*3: $MAX_DEPTH" 
#APPLY MEAN DEPTH AND MAX DEPTH FILTERS
#ALSO in ONE GO APPLY THE Filter 2: MASKING SPURIOUS HETEROZYGOTES
#  filter imbalanced het genotypes with binomial test
#     p-value of 0.01 tighten to 0.0001 for more conservative (keep more borderline hets)
#     for example if a genotype has meanDepth ~34 reads and has 5 reads of minorAllele will fail 
#     or loosen to 0.05 if you want to be stricter about balance.
#NOTE: Don't rely on GT masking to protect your ANGSD/NGSadmix results, will still read GL PL, but masking GT is still good for vcf based analysis and missingness filters
bcftools view -e "INFO/DP / INFO/NS > ${MAX_DEPTH} || INFO/DP / INFO/NS < 3" ${prefix}.temp1.vcf.gz | \
bcftools +setGT -- -t "b:AD<0.0001" -n . | \
bgzip -c > ${prefix}.temp2.vcf.gz
bcftools index ${prefix}.temp2.vcf.gz
#FILTER PARALOGS (LOCI WITH TOO MANY SNPS)
#count number of SNPs per loci
bcftools query -f '%CHROM\n' ${prefix}.temp2.vcf.gz \
| awk '{count[$1]++} END {for (locus in count) print count[locus], locus}' \
| sort -k1,1n \
> ${prefix}.temp2.numberSNPs_per_loci.txt
awk '{dist[$1]++} END {for (n in dist) print dist[n], n}' ${prefix}.temp2.numberSNPs_per_loci.txt \
| sort -k2,2n \
> ${prefix}_numberSNPs_per_loci_distribution.txt
#For rad seq > 4 (population level) is suspicious, but in this genus level assembly we many more.
#Distribution is not bimodal —> no clean trough separating "real loci" from "paralogs". 
#whitelist loci with less than 10 SNPs (the rest are likely paralog piled up) after which we see sharp decline
awk '$1 <= 8 {print $2}' ${prefix}.temp2.numberSNPs_per_loci.txt > ${prefix}.temp2.keep_loci.txt
#Since each contig spans positions 1–36, create a proper regions file:
awk '{print $1"\t1\t36"}' ${prefix}.temp2.keep_loci.txt > ${prefix}.temp2.regions.txt
#subset vcf by "contigs" (radloci)
bcftools view -R ${prefix}.temp2.regions.txt ${prefix}.temp2.vcf.gz | bgzip -c > ${prefix}.vcf.gz
bcftools index ${prefix}.vcf.gz
echo "check how many variants remain"
bcftools query -f '%CHROM' ${prefix}.vcf.gz | wc -l
echo "check how many loci remain"
bcftools query -f '%CHROM' ${prefix}.vcf.gz | sort | uniq | wc -l
#compute conprehensive stats
bcftools stats -s - ${prefix}.vcf.gz > ${prefix}.stats

#EOF SCRIPT
#PUT ON A ASCRIPT AND RUN
# FIRST Prepare input (compress + index)
# bgzip populations_n4.snps.vcf
# bcftools index populations_n4.snps.vcf.gz
bsub -G compute-christine.e.edwards -q general -M 16GB -R "rusage[mem=16GB]" -a "docker(staphb/bcftools)" "./bcftools_filters.sh populations_n3.snps.vcf.gz all_n3_wDP_AD0001"
bsub -G compute-christine.e.edwards -q general -M 16GB -R "rusage[mem=16GB]" -a "docker(staphb/bcftools)" "./bcftools_filters.sh populations_n4.snps.vcf.gz all_n4_wDP_AD0001"
bsub -G compute-christine.e.edwards -q general -M 16GB -R "rusage[mem=16GB]" -a "docker(staphb/bcftools)" "./bcftools_filters.sh populations_n5.snps.vcf.gz all_n5_wDP_AD0001"
#download from cluster all assembly outputs
#scp -r a.eduardo@compute1-client-1.ris.wustl.edu:/storage1/fs1/christine.e.edwards/Active/weinmannia/assembly/populations .
