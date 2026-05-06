#!/bin/bash
INPUT=$1
OUTPUT_PATH=$2
seqkit grep -P --degenerate -R 13:24 -s -p 'CGANNNNNNTGC' $INPUT -o temp1${INPUT}
seqkit grep -P --degenerate -R 13:24 -s -p 'GCANNNNNNTCG' $INPUT | seqkit seq --complement --reverse -o temp2${INPUT}
cat temp1${INPUT} temp2${INPUT} > $OUTPUT_PATH/$INPUT


#### PUT THE ABOVE LINES INTO A SCRIPT AND MAKE IT EXECUTABLE THEN RUN THE FOLLOWINNG COMMAND ON A DIRECTORY WITH SAMPLES
#!/bin/bash
for fastq in $(ls *fastq.gz)
do
./normalize_strands.sh $fastq /scratch1/fs1/christine.e.edwards/eduardo_tmp/samples/normalized
done


export SCRATCH1=/scratch1/fs1/christine.e.edwards
export STORAGE1=/storage1/fs1/christine.e.edwards/Active
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"
bsub -G compute-christine.e.edwards -q general -a "docker(staphb/seqkit)" -R "rusage[mem=16GB]" /scratch1/fs1/christine.e.edwards/eduardo_tmp/samples/run_orient.sh



# for seq in $(ls | grep -v temp | grep fastq.gz)
# do
# cat temp1$seq temp2$seq > ./normalized/$seq
# done