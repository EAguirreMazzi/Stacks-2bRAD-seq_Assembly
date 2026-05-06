#!/bin/bash
cd /scratch1/fs1/christine.e.edwards/eduardo_tmp/samples
FASTP="/storage1/fs1/christine.e.edwards/Active/aguirre/fastp"
NUM_JOBS=10

files=$(ls | grep ".fastq")
rm fastp_cmds
for file in $files
do
  base=$(basename $file .fastq)
  echo "$FASTP -i ${file} -o ${file}.gz -q 20 -u 20 --length_required 36 --trim_poly_g --poly_g_min_len 6 --n_base_limit 0 -h ${base}.html -j ${base}.json" >> fastp_cmds
done

split -d -n l/$NUM_JOBS fastp_cmds fastp_cmd_
#SUMBIT JOBS

export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"

for fastp_cmd_ in $(ls | grep fastp_cmd_); do
  sed -i '1i #!/bin/bash' $fastp_cmd_
  chmod +x $fastp_cmd_
  bsub -G compute-christine.e.edwards -g /a.eduardo/eaguirre -q general -a "docker(continuumio/miniconda3)" -R "rusage[mem=8GB]" ./$fastp_cmd_
done
