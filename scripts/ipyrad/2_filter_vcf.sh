bgzip phylo_base.vcf
bcftools index phylo_base.vcf.gz

input_vcf=phylo_base.vcf.gz
prefix=phylo_base_clean
#FILTER PARALOGS (LOCI WITH TOO MANY SNPS)
#count number of SNPs per loci
# bcftools query -f '%CHROM\n' ${input_vcf} \
# | awk '{count[$1]++} END {for (locus in count) print count[locus], locus}' \
# | sort -k1,1n \
# > ${prefix}.temp2.numberSNPs_per_loci.txt
# awk '{dist[$1]++} END {for (n in dist) print dist[n], n}' ${prefix}.temp2.numberSNPs_per_loci.txt \
# | sort -k2,2n \
# > ${prefix}_numberSNPs_per_loci_distribution.txt
# #For rad seq > 4 (population level) is suspicious, but in this genus level assembly we many more.
# #Distribution is not bimodal —> no clean trough separating "real loci" from "paralogs". 
# #whitelist loci with less than 10 SNPs (the rest are likely paralog piled up) after which we see sharp decline
# awk '$1 <= 8 {print $2}' ${prefix}.temp2.numberSNPs_per_loci.txt > ${prefix}.temp2.keep_loci.txt
# #Since each contig spans positions 1–36, create a proper regions file:
# awk '{print $1"\t1\t36"}' ${prefix}.temp2.keep_loci.txt > ${prefix}.temp2.regions.txt
# #subset vcf by "contigs" (radloci)
# bcftools view -R ${prefix}.temp2.regions.txt ${prefix}.temp2.vcf.gz | bgzip -c > ${prefix}.vcf.gz
# bcftools index ${prefix}.vcf.gz


MIN_MAC=1

for MAX_F_MISSING in 0.8 0.7 0.6 0.5
do
prefix=phylo_base_FMIS${MAX_F_MISSING}
bcftools view -i "F_MISSING<=${MAX_F_MISSING}" $input_vcf \
  | bcftools view -i "MAC>${MIN_MAC}" \
  | bgzip -c > ${prefix}.vcf.gz
bcftools index ${prefix}.vcf.gz
done

vcf2phylip=/storage1/fs1/christine.e.edwards/Active/aguirre/vcf2phylip/vcf2phylip.py
for MAX_F_MISSING in 1 0.9
do
prefix=phylo_base_FMIS${MAX_F_MISSING}
  $vcf2phylip -i ${prefix}.vcf.gz
  mv ${prefix}.min4.phy ${prefix}.phy
done


$vcf2phylip -i phylo_base.vcf.gz -m 10
$vcf2phylip -i phylo_base.vcf.gz -m 15
$vcf2phylip -i phylo_base.vcf.gz -m 20
$vcf2phylip -i phylo_base.vcf.gz -m 25
$vcf2phylip -i phylo_base.vcf.gz -m 30
$vcf2phylip -i phylo_base.vcf.gz -m 35
$vcf2phylip -i phylo_base.vcf.gz -m 40
$vcf2phylip -i phylo_base.vcf.gz -m 45
$vcf2phylip -i phylo_base.vcf.gz -m 50



for msl in 10 15 20 25 30 35 40 45 50; do
iqtree -s phylo_base.min${msl}.phy -m GTR+ASC -T 8 --prefix phylo_base.min${msl}
done



for msl in 10 15 20 25 30 35 40 45 50; do
bsub -G compute-christine.e.edwards -g /a.eduardo/eaguirre -q general -R "rusage[mem=8GB]" -a "docker(condaforge/mambaforge)" bash \
        -c "conda init bash && source $HOME/.bashrc && conda activate ipyrad_0.9.90 && raxmlHPC -f d -m ASC_GTRCAT --asc-corr=lewis -n phylo_base.min${msl} -s phylo_base.min${msl}.varsites.phy -p 12345"
done

for msl in 10 15 20 25 30 35 40 45 50; do
bsub -G compute-christine.e.edwards -g /a.eduardo/eaguirre -q general -R "rusage[mem=8GB]" -a "docker(condaforge/mambaforge)" bash \
        -c "conda init bash && source $HOME/.bashrc && conda activate iqtree && iqtree -s phylo_base.min${msl}.varsites.phy -m GTR+ASC -T 4 --prefix phylo_base.min${msl}"
done









bkill 858610
bkill 858611
bkill 858612
bkill 858613
bkill 858614
bkill 858615
bkill 858616
bkill 858617


raxmlHPC -d -s phylo_base.min50.varsites.phy -m ASC_GTRCAT --asc-corr=lewis -n phylo_base.min50 -p 12345