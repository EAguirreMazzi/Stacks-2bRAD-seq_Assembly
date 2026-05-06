#COPY ALL RAW FILES TO ONE DIRECTORY AND SIMPLYFY NAME FOR EASIER MANIPULATION LATER
#Easier For plates 1 to 14 no plate index (one plate per plate)
cd /storage1/fs1/christine.e.edwards/Active/weinmannia/data
mkdir -p all_raw
for d in Plate{1..14}; do
  [ -d "$d" ] || continue

  for f in "$d"/raw_data/*[Rr][Oo][Ww]*.fastq.gz; do
    [ -e "$f" ] || continue

    base=$(basename "$f")

    # extract row number (case-insensitive)
    if [[ $base =~ [Rr][Oo][Ww]([0-9]+) ]]; then
      row="${BASH_REMATCH[1]}"
      new="all_raw/${d}_row${row}.fastq.gz"

      cp -n -- "$f" "$new"
    fi
  done
done
#for Plates 15 and 16 need to add the PCR barcode label
#Ill do manually since BC labels are not consistent and I'll use only R1 (this is very important to avoid batch effects downstream)
#For Plate15 only TS2 and TS3 have weinmannia data so will use only this ones
cd /storage1/fs1/christine.e.edwards/Active/weinmannia/data/Plate15/raw_data
ALL_RAW=/storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw
cp LeiWeiAri_PL2_ROW1_S9_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row1.fastq.gz
cp LeiWeiAri_PL2_ROW2_S10_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row2.fastq.gz
cp LeiWeiAri_PL3_ROW1_S17_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row1.fastq.gz
cp LeiWeiAri_PL3_ROW2_S18_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row2.fastq.gz
cp LeiWeiAri_PL2_ROW3_S11_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row3.fastq.gz
cp LeiWeiAri_PL2_ROW4_S12_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row4.fastq.gz
cp LeiWeiAri_PL3_ROW3_S19_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row3.fastq.gz
cp LeiWeiAri_PL3_ROW4_S20_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row4.fastq.gz
cp LeiWeiAri_PL2_ROW5_S13_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row5.fastq.gz
cp LeiWeiAri_PL2_ROW6_S14_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row6.fastq.gz
cp LeiWeiAri_PL3_ROW5_S21_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row5.fastq.gz
cp LeiWeiAri_PL3_ROW6_S22_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row6.fastq.gz
cp LeiWeiAri_PL2_ROW7_S15_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row7.fastq.gz
cp LeiWeiAri_PL2_ROW8_S16_R1_001.fastq.gz $ALL_RAW/Plate15_PL2_row8.fastq.gz
cp LeiWeiAri_PL3_ROW7_S23_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row7.fastq.gz
cp LeiWeiAri_PL3_ROW8_S24_R1_001.fastq.gz $ALL_RAW/Plate15_PL3_row8.fastq.gz
# Same for PLATE 16 Labeling and filesystem is not consistent with Plate 15 so will do manually
# also will merge two runs one to compensate subsampling of loci
ALL_RAW=/storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw
RUN1=/storage1/fs1/christine.e.edwards/Active/weinmannia/data/Plate16/run_1_with_batch_effect/01.RawData
RUN2=/storage1/fs1/christine.e.edwards/Active/weinmannia/data/Plate16/01.RawData
cat $RUN1/Wein_16_row01/Wein_16_row01_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row1/Po2We16_TS1_row1_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row1.fastq.gz
cat $RUN1/Wein_16_row02/Wein_16_row02_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row2/Po2We16_TS1_row2_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row2.fastq.gz
cat $RUN1/Wein_16_row03/Wein_16_row03_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row3/Po2We16_TS1_row3_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row3.fastq.gz
cat $RUN1/Wein_16_row04/Wein_16_row04_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row4/Po2We16_TS1_row4_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row4.fastq.gz
cat $RUN1/Wein_16_row05/Wein_16_row05_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row5/Po2We16_TS1_row5_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row5.fastq.gz
cat $RUN1/Wein_16_row06/Wein_16_row06_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row6/Po2We16_TS1_row6_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row6.fastq.gz
cat $RUN1/Wein_16_row07/Wein_16_row07_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row7/Po2We16_TS1_row7_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row7.fastq.gz
cat $RUN1/Wein_16_row08/Wein_16_row08_CKDL260001793-1A_23575FLT4_L1_1.fq.gz $RUN2/Po2We16_TS1_row8/Po2We16_TS1_row8_CKDL260006557-1A_23737JLT4_L4_1.fq.gz > $ALL_RAW/Plate16_TS1_row8.fastq.gz


#THEN MOVE TO SCRATCH FOR BETTER PERFORMANCE
SCRATCH_DIR=/scratch1/fs1/christine.e.edwards/eduardo_tmp/all_raw
mkdir -p SCRATCH_DIR
rsync -a --info=progress2 /storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw/* $SCRATCH_DIR