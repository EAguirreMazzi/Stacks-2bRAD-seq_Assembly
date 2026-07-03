
#get the list (from n3 admixture analysis)
#choose n3 beacuase more stable phylogenetic placement of taxa
#also reasonable parameter accorfding to segovia et al. 2025 ct0.92 (which we will use here too)
bcftools query -l ../3_evolution/sequence_data/phylogenomics/phylo_v1/n3_phyV1_FMIS0.5_MAC0.vcf.gz > scripts/ipyrad/phylo_base_sample_list.txt
#
mkdir -p /scratch1/fs1/christine.e.edwards/eduardo_tmp/phylo_base/fastq_files/
for seqID in $(cat scripts/ipyrad/phylo_base_sample_list.txt); do
echo "rsync -a /storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw/normalized/${seqID}.fastq.gz /scratch1/fs1/christine.e.edwards/eduardo_tmp/phylo_base/fastq_files/"
done > scripts/ipyrad/move_to_scratch.txt
#now samples are in scratch lets start the pipeline step by step to control for errors and pre-filter before clustering

cd /scratch1/fs1/christine.e.edwards/eduardo_tmp/phylo_base
#create new assembly
ipyrad -n phylo_base
## edit the params file to enter your raw_fastq_path and barcodes path
sed -i '/\[4] /c\/scratch1/fs1/christine.e.edwards/eduardo_tmp/phylo_base/fastq_files/*fastq.gz  ## [4] ' params-phylo_base.txt
sed -i '/\[7] /c\2brad  ## [7] ' params-phylo_base.txt
#depth parametes 
sed -i '/\[11] /c\5  ## [1] ' params-phylo_base.txt
sed -i '/\[12] /c\3  ## [12] ' params-phylo_base.txt
sed -i '/\[13] /c\2000  ## [13] ' params-phylo_base.txt
#clustering threshold (withins samples step3)
sed -i '/\[14] /c\0.9444444444444444  ## [14] ' params-phylo_base.txt
sed -i '/\[17] /c\36  ## [17] ' params-phylo_base.txt
sed -i '/\[19] /c\0  ## [19] ' params-phylo_base.txt
sed -i '/\[20] /c\0.1  ## [20] ' params-phylo_base.txt
sed -i '/\[22] /c\0.2  ## [22] ' params-phylo_base.txt
sed -i '/\[23] /c\0  ## [23] ' params-phylo_base.txt
sed -i '/\[27] /c\*  ## [27] ' params-phylo_base.txt

#RUN for all
ipyrad -p params-phylo_base.txt -s12 -c8
#then subset batches for efficiency (will be merged before S6)
ipyrad -p params-phylo_base.txt -b phyloV2_1 agui_49_1 fuen_21911_8 fuen_19134_1 fuen_19973_1 fuen_22286_1 fuen_22037_1 cayo_4867_1 fuen_22125_1 chas_303_1 chas_339_1 chas_523_1 loza_1786_1 fuen_22249_1 fuen_21995_1 fuen_22204_1 fuen_22224_1 chas_508_1 chas_239_1 chas_279_1 vasq_46388_1 cast_2140_8
ipyrad -p params-phylo_base.txt -b phyloV2_2 chas_244_1 fuen_21862_1 corn_1181_1 fuen_21706_1 fuen_22290_1 fuen_22017_1 mald_3229_1 sall_3086_1 chas_301_1 chas_447_1 fuen_22152_1 loza_1617_1 fuen_22275A_1 vasq_46389_1 fuen_22206_1 fuen_22227_8 chas_458_1 fuen_20089_1 chas_280_1 agui_50_1 idar_8002_1
ipyrad -p params-phylo_base.txt -b phyloV2_3 mcph_15916_1 fuen_21856_1 fuen_19205_1 fuen_21705_1 fuen_22342_1 fuen_22039_1 fuen_22099_1 fuen_22356_1 chas_536_1 chas_347_1 fuen_21918_8 fuen_22352_1 fuen_22287_1 cayo_4428_1 fuen_22183_1 fuen_22066_8 chas_468_1 chas_453_1 chas_248_1 agui_45_1 rhum_158_1
ipyrad -p params-phylo_base.txt -b phyloV2_4 agui_2_1 fuen_21857_1 chas_615_8 fuen_21698_1 fuen_22344_1 fuen_22049_1 fuen_20098_1 fuen_22077_1 faja_4626_1 neil_15248_1 fuen_21912_1 cach_91_1 cabe_P225_1 fuen_22331A_1 fuen_22200_1 fuen_19421_1 chas_457_1 calb_297_8 chas_673_8 agui_15_1 rhum_157_1
ipyrad -p params-phylo_base.txt -b phyloV2_5 agui_1_1 fuen_21907_1 fuen_21775_1 fuen_21699_1 fuen_22347_1 fuen_10090_1 fuen_20093_1 faja_6376_1 chas_623_1 chas_220_1 fuen_21910_1 sall_2943_1 fuen_22255_1 fuen_22294_1 fuen_22243_1 fuen_22032_1 fuen_21949_1 chas_403_1 mash_134B_1 mcph_18115_1 idar_7973_1
ipyrad -p params-phylo_base.txt -b phyloV2_6 agui_75_1 fuen_19149_1 fuen_21932_1 fuen_22312_1 sall_2992_1 fuen_22038_1 fuen_20085_1 agui_52_1 chas_389_1 fuen_22233_1 fuen_22317A_1 fuen_22317B_1 sall_2850_8 fuen_22297A_1 chas_1003_1 fuen_22044_1 fuen_20051_1 chas_384_8 mash_220_1 lowr_5713_1 
ipyrad -p params-phylo_base.txt -b phyloV2_7 agui_32_1 fuen_16131_1 fuen_21788_1 fuen_21715_1 fuen_22251_1 fuen_10442_1 fuen_22167_1 chas_431_1 fuen_20064A_1 agui_35_1 fuen_21885_1 fuen_22298_1 fuen_19892_1 fuen_21720_1 chas_263_1 fuen_22176_1 fuen_22348_8 chas_663_1 chas_666_1 raso_125_1 
ipyrad -p params-phylo_base.txt -b phyloV2_8 chas_234_2 fuen_19138_1 fuen_21940_1 fuen_22316A_1 sall_3020_1 fuen_14013_1 fuen_21714_1 chas_501_1 fuen_20050_1 agui_30_1 fuen_21868_1 fuen_22262_1 fuen_21710_8 corn_1125_1 chas_269_8 fuen_22160_8 sall_2921_1 chas_411_1 chas_277_1 rako_179_1 
ipyrad -p params-phylo_base.txt -b phyloV2_9 agui_33_1 fuen_19165_1 cayo_3710_1 fuen_20053_1 vale_8845_1 fuen_13947_1 fuen_20059_1 chas_419_1 fuen_22016_1 agui_72_1 fuen_21865_1 fuen_22250_1 fuen_21708_8 fuen_21780_1 fuen_22221_1 sall_2895_8 chas_428_1 sall_2923_1 chas_474_1 orti_607_1 
ipyrad -p params-phylo_base.txt -b phyloV2_10 agui_20_1 fuen_19143_1 agui_69_1 fuen_22315B_1 fuen_22328_1 maju_8773_1 fuen_20095_1 chas_511_8 sall_3059_8 agui_74_1 sego_4_1 fuen_22292_1 fuen_21717_8 sall_3032_1 fuen_22228_1 faja_6365_8 chas_482_1 chas_401_1 chas_580_1 orti_602_1 
ipyrad -p params-phylo_base.txt -b phyloV2_11 agui_54_1 fuen_21982_1 sall_3045_1 fuen_22303_1 fuen_22238_1 fuen_22085_1 vasq_46393_1 chas_510_1 sall_2898_1 agui_56_1 sego_1_8 fuen_22267_1 fuen_21686_8 sall_3024_1 fuen_22013_8 sall_2892_1 chas_541_1 chas_408_1 cabe_P228_8 faja_6418_8
ipyrad -p params-phylo_base.txt -b phyloV2_12 chas_549_1 fuen_21872_1 chas_509_1 fuen_22324_1 fuen_22201_1 fuen_22150_1 fuen_20069_1 chas_274_1 chas_331_1 agui_77_1 cayo_4887_1 vill_106_1 fuen_22258_1 sall_3027_1 fuen_22215_1 fuen_22073_1 chas_621_8 chas_449_1 sall_2913_1 lueb_3775_8 
ipyrad -p params-phylo_base.txt -b phyloV2_13 fuen_19215_1 neil_16242_1 cayo_3779_1 corn_1102_1 fuen_22218_1 fuen_20060_1 fuen_22164_8 chas_427_1 chas_337_1 agui_81_1 fuen_22003_1 fuen_22330_1 fuen_22299_1 sall_3061_8 fuen_22226_1 fuen_20097_1 sall_2851_1 chas_264_1 chas_343_1 calb_299_1 
ipyrad -p params-phylo_base.txt -b phyloV2_14 chas_635_8 fuen_22143_1 fuen_19274_1 fuen_21978_1 fuen_22282_1 fuen_22223_1 chas_472_1 faja_6377_1 chas_349_1 agui_80_1 fuen_22001_1 fuen_19239_1 fuen_19858_1 fuen_19923_1 sall_2901_1 fuen_22354_1 chas_445_1 calb_298_1 chas_344_1 faja_5941_1 
ipyrad -p params-phylo_base.txt -b phyloV2_15 fuen_21915_1 fuen_21914_1 sall_3036_8 fuen_22304_1 sall_3002_1 fuen_22154C_1 fuen_22173_8 chas_421_1 chas_593_1 agui_36_1 sego_3_8 fuen_22261_1 fuen_21740_1 sall_3030_8 fuen_19115_1 chas_371_1 sall_2908_2 chas_405_8 chas_693_1 faja_6423_8 

#always run this block before sending jobs within scratch
export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"
export CONDA_PKGS_DIRS=/storage1/fs1/christine.e.edwards/Active/aguirre/conda/pkgs/
export CONDA_ENVS_DIRS=/storage1/fs1/christine.e.edwards/Active/aguirre/conda/envs/
# RUn steps 3 to 6 for each batch in parallel
for batch in 11; do 
  bsub \
    -G compute-christine.e.edwards \
    -g /a.eduardo/eaguirre \
    -q general \
    -n 4 \
    -R "rusage[mem=32GB]" \
    -J ipy${batch} \
    -a "docker(condaforge/mambaforge)" \
    bash -c "conda init bash && source $HOME/.bashrc && conda activate ipyrad_0.9.90 && ipyrad -p params-phyloV2_${batch}.txt -s 345 -c 4 -t 1 -f"
done


ipyrad -p params-phyloV2_15.txt -r

#merge all in one then run step 6
ipyrad -m phyloV2 params-phyloV2_1.txt params-phyloV2_2.txt params-phyloV2_3.txt params-phyloV2_4.txt params-phyloV2_5.txt params-phyloV2_6.txt params-phyloV2_7.txt params-phyloV2_8.txt params-phyloV2_9.txt params-phyloV2_10.txt params-phyloV2_11.txt params-phyloV2_12.txt params-phyloV2_13.txt params-phyloV2_14.txt params-phyloV2_15.txt -f

#change clustering threshold across samples (genus level expectation)
sed -i '/\[14] /c\0.8888889  ## [14] ' params-phyloV2.txt


# run the job
export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"
export CONDA_PKGS_DIRS=/storage1/fs1/christine.e.edwards/Active/aguirre/conda/pkgs/
export CONDA_ENVS_DIRS=/storage1/fs1/christine.e.edwards/Active/aguirre/conda/envs/


bsub \
    -G compute-christine.e.edwards \
    -g /a.eduardo/eaguirre \
    -q general \
    -n 24 \
    -R "rusage[mem=128GB]" \
    -J ipyrads67 \
    -a "docker(condaforge/mambaforge)" \
    bash -c "conda init bash && source $HOME/.bashrc && conda activate ipyrad_0.9.90 && ipyrad -p params-phyloV2.txt -s 67 -c 24 -t 1 -f"