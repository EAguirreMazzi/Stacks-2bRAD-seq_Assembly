#!/bin/bash
cd /scratch1/fs1/christine.e.edwards/eduardo_tmp/all_raw
TRIM="/storage1/fs1/christine.e.edwards/Active//aguirre/scripts_preanalisis/trim2bRAD_2barcodes_noAdap.pl"
NUM_JOBS=20
files=$(ls | grep ".fastq" | grep -v "done_raws")
rm trim_cmds
for file in $files
do
  #echo "perl ${TRIM} fastq=$file site='.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}' barcode2='[ATGC]{4}'" >> trim_cmds
  #filter sites with no NG adaptor
  echo "perl ${TRIM} fastq=$file site='.{1}G.{10}CGA.{6}TGC.{10}C.{1}|.{1}G.{10}GCA.{6}TCG.{10}C.{1}' barcode2='[ATGC]{4}'" >> trim_cmds
done
split -d -n l/$NUM_JOBS trim_cmds trim_cmd_
rm trim_cmds
#SUMBIT JOBS
export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"

for trim_cmd_ in $(ls | grep trim_cmd_); do
  sed -i '1i #!/bin/bash' $trim_cmd_
  chmod +x $trim_cmd_
  #bsub -G compute-christine.e.edwards -g /a.eduardo/eaguirre -q general -a "docker(condaforge/mambaforge)" -R "rusage[mem=8GB]" ./$trim_cmd_
done
