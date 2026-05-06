#!/bin/bash
cd /scratch1/fs1/christine.e.edwards/eduardo_tmp/samples/
USTACKS="/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/ustacks"
NUM_JOBS=10

mkdir -p /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly

files=$(ls | grep ".fastq")

for file in $files
do
  base=$(basename $file .fastq)
  echo "$USTACKS -f ${file} -o /scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly -m 3 -M 2 --model-type bounded --bound-low 0.001 --bound-high 0.07 --alpha 0.05 --disable-gapped -t 8" >> ustacks_cmds
done

split -d -n l/$NUM_JOBS ustacks_cmds ustacks_cmd_
#SUMBIT JOBS

export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"

for ustacks_cmd_ in $(ls | grep ustacks_cmd_); do
  #sed -i '1i #!/bin/bash' $ustacks_cmd_
  #chmod +x $ustacks_cmd_
  bsub -o ${ustacks_cmd_}.log -G compute-christine.e.edwards -g /a.eduardo/eaguirre -q general -a "docker(gcc:12.2.0)" -R "rusage[mem=16GB]" ./$ustacks_cmd_
done


#After jobs finished extract information from (stdout)log files.
cat ustacks_cmd_*.log > ustacks_stdout.log
logfile=ustacks_stdout.log
header="Sequence_ID\tloci_stacks\tmean_coverage\tsd_coverage\tmax_coverage\tn_reads\tpercentage_reads"
paste <(cat $logfile | grep "Input file:") <(cat $logfile | grep "Final number of stacks built: ") <(cat $logfile | grep "Final coverage: ") > ustacks_coverage_stats.txt
sed -i "s/  Input file: '//" ustacks_coverage_stats.txt
sed -i "s/.fastq.gz'\tFinal number of stacks built: /\t/" ustacks_coverage_stats.txt
sed -i "s/\tFinal coverage: mean=/\t/" ustacks_coverage_stats.txt
sed -i "s/; stdev=/\t/" ustacks_coverage_stats.txt
sed -i "s/; max=/\t/" ustacks_coverage_stats.txt
sed -i "s/; n_reads=/\t/" ustacks_coverage_stats.txt
sed -i "s/(/\t/" ustacks_coverage_stats.txt
sed -i "s/%)//" ustacks_coverage_stats.txt


cat $header ustacks_coverage_stats.txt > ustacks_coverage_stats.txt
