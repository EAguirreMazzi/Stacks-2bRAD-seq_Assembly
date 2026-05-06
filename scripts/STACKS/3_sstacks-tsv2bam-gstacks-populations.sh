#!/bin/bash
#BSUB -G compute-christine.e.edwards
#BSUB -q general
#BSUB -n 16
#BSUB -R "rusage[mem=128GB]"
#BSUB -a "docker(gcc:12.2.0)"
# Limit hidden threading (important)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
SCRATCH_DIR=/scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly
POPMAP=/scratch1/fs1/christine.e.edwards/eduardo_tmp/assembly/raw_full_popmap.txt
POP_OUT=$SCRATCH_DIR/populations/n3/raw_full
/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/sstacks -P $SCRATCH_DIR -M $POPMAP --disable-gapped -t 16
/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/tsv2bam -P $SCRATCH_DIR -M $POPMAP -t 16
/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/gstacks -P $SCRATCH_DIR -M $POPMAP --kmer-length 21 --max-debruijn-reads 500 --min-kmer-cov 1 -t 16
mkdir -p $POP_OUT
/storage1/fs1/christine.e.edwards/Active/aguirre/stacks/stacks-2.68/populations -P $SCRATCH_DIR -M $POPMAP -O $POP_OUT --no-hap-exports --vcf-all -t 16