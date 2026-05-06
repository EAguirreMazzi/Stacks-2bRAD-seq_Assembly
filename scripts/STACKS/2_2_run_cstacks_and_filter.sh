#!/bin/bash
max_missmatches=$1
OUTPATH=/scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog_n${max_missmatches}
mkdir -p $OUTPATH
/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/cstacks -n ${max_missmatches} \
    -s ./chas_510_1 -s ./chas_356_1 -s ./idar_8002_1 -s ./chas_522_1 -s ./cayo_4939_1 -s ./agui_1_1 \
    -s ./mcph_18115_1 -s ./fuen_19858_1 -s ./chas_445_1 -s ./cayo_4940_1 -s ./maju_8773_1 -s ./fuen_19441_1 \
    -s ./fuen_22164_8 -s ./chas_666_1 -s ./sall_3030_8 -s ./faja_6400_8 -s ./sall_2923_1 -s ./agui_80_1 \
    -s ./chas_419_1 -s ./fuen_19165_1 -s ./cayo_4909_1 -s ./fuen_21914_1 -s ./cabe_P225_1 -s ./chas_250_1 \
    -s ./idar_7982_8 -s ./faja_6391_1 -s ./brad_1089_1 -s ./sall_3044_8 -s ./chas_689_1 -s ./fuen_19205_1 \
    -s ./lueb_3775_8 -s ./fuen_22003_1 -s ./rako_179_1 -s ./fuen_21731_1 -s ./chas_457_1 -s ./fuen_22099_1 \
    -s ./fuen_22344_1 -s ./fuen_19240_1 -s ./sall_3033_8 -s ./mcph_15916_1 -s ./cayo_4930_1 -s ./fuen_19204_1 \
    -s ./sall_2898_1 -s ./fuen_22311_1 -s ./sall_2996_1 -s ./sall_2858_1 -s ./vale_8845_1 -s ./fuen_21687_1 \
    -s ./samo_11_1 -s ./chas_528_8 -s ./rhum_158_1 -s ./chas_411_1 -s ./sall_3058_8 -s ./agui_19_1 -s ./fuen_19143_1 \
    -s ./neil_15248_1 -s ./fuen_21712_1 -s ./fuen_22001_1 -s ./chas_462_1 -s ./fuen_21885_1 -s ./fuen_20091_1 \
    -s ./chas_549_1 -s ./chas_523_1 -s ./faja_6365_8 -s ./chas_272_1 -s ./fuen_22304_1 -s ./fuen_22348_8 -s ./chas_438_1 \
    -s ./chas_283_1 -s ./calb_296_1 -s ./fuen_22035_1 -s ./chas_562_1 -s ./vasq_46389_1 -s ./fuen_21724_1 -s ./chas_555_1 \
    -s ./fuen_19973_1 -s ./faja_6389_8 -s ./lowr_5713_1 -s ./chas_505_1 -s ./agui_50_1 -s ./sego_4_1 -s ./chas_468_1 \
    -s ./faja_6423_8 -s ./sego_1_8 -s ./fuen_22053_1 -s ./cayo_4867_1 -s ./orti_602_1 -s ./sall_2864_1 -s ./chas_624_1 \
    -s ./cayo_4935_1 -s ./chas_453_1 -s ./fuen_15146_1 -s ./fuen_21656_1 -s ./agui_81_1 \
    -o $OUTPATH --disable-gapped -t 12



#COPY THE ABOVE BLOCK to A SCRIPT THEN RUN:

# Limit hidden threading (important), more explicit use of threading within stacks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"

cd /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly
for n in {3..6}
do
OUTPATH=/scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog_n${n}
mkdir -p $OUTPATH
bsub -o ${OUTPATH}/catalo_n${n}.log -G compute-christine.e.edwards -q general -R "rusage[mem=32GB]" -a "docker(gcc:12.2.0)" "./cstacks_cmd.sh ${n}"
done

#or run with popmaps not good for branching but convenient
#/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/cstacks \
#    -P /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly \
#    -n 3 -M /scratch1/fs1/christine.e.edwards/eduardo_tmp/catalog_builders_popmap.txt --disable-gapped -t 12
#cp /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog.* /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog_n3

#for i in {2..3}
#do
#/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/cstacks -P /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly --catalog /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog -n 3 -M /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/popmaps/catalog_builders_stage${i}_popmap.txt --disable-gapped -t 12
#cp /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog.* /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/catalog_n3
#echo "stage ${i} finished and backed up" >> cstacks_progress.log
#done

# EXPLORE SAMPLE OCURRENCE IN CATALOG LOCI
# For each loci get
#1. Total number of matches (list length)
#2. Number of UNIQUE samples matching (unique values)
#3. Number of samples with MULTIPLE matches (repeated values)
#4. Maximum stacks per sample at this locus (ex. 5,5,5,5,3,3,1 = 4)

zcat catalog.tags.tsv.gz | grep -v "#" | awk -F'\t' '
{  
    # If locus ID is in column 1:
    locus_id = $1
    
    # Process column 4 (tags)
    split($4, tags, ",")
    
    # Count samples for this locus
    delete count
    for(i in tags) {
        split(tags[i], read_parts, "_")
        sample = read_parts[1]
        count[sample]++
    }
    
    # Calculate metrics
    total = 0
    unique = 0
    multiple = 0
    max_cnt = 0
    
    for(s in count) {
        total += count[s]
        unique++
        if(count[s] > 1) multiple++
        if(count[s] > max_cnt) max_cnt = count[s]
    }
    
    # Output locus ID and metrics
    printf "%s\t%d\t%d\t%d\t%d\n", locus_id, total, unique, multiple, max_cnt
}' > catalog.occurrence.txt
#ecplore distributions of this metrics to decide particularly the max_cnt
cut -f 5 catalog.occurrence.txt | sort -n | uniq -c
#whitelist loci to keep 
#         .max_cnt   .min_samples
awk '$5 < 5 && $3 >= 3 {print $1}' catalog.occurrence.txt > whitelist_loci.txt
mkdir -p ./filtered_catalog
zcat catalog.tags.tsv.gz | awk 'NR==FNR {whitelist[$1]; next} ($1 in whitelist)' whitelist_loci.txt - | gzip > ./filtered_catalog/catalog.tags.tsv.gz
zcat catalog.snps.tsv.gz | awk 'NR==FNR {whitelist[$1]; next} ($1 in whitelist)' whitelist_loci.txt - | gzip > ./filtered_catalog/catalog.snps.tsv.gz
zcat catalog.alleles.tsv.gz | awk 'NR==FNR {whitelist[$1]; next} ($1 in whitelist)' whitelist_loci.txt - | gzip > ./filtered_catalog/catalog.alleles.tsv.gz


#explore distribution of filtered catalog to make sure it matches expected
cd filtered_catalog
cut -f 5 catalog.occurrence.txt | sort -n | uniq -c
#also figure out if there is not Batch effect
zcat catalog.tags.tsv.gz | cut -f 4 > catalog.tags.column4.txt
sed -i 's/_[0-9]*//g' catalog.tags.column4.txt
sed 's/,/\n/g' catalog.tags.column4.txt | sort -n | uniq -c > loci_per_sample.txt
#explore this table (make more explicit next time so ID is hardcoded matched)
paste <(cat loci_per_sample.txt) <(zcat ../catalog.sample_list.tsv.gz | grep -v "#")