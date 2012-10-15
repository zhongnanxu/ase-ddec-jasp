#!/bin/sh
###################################################
#PBS -N matlab 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00
#PBS -l pmem=3900mb
##PBS -l vmem=1gb
#PBS -q joe-medium
#PBS -r n
#PBS -j oe
#PBS -V
###################################################

cd $PBS_O_WORKDIR

echo "Nodes chosen are:"
cat $PBS_NODEFILE

EXE=/usr/local/packages/matlab/bin/matlab

$EXE -nodesktop -r calculate_nonperiodic_RMAE > nonperiodic_RMAE_output.txt
