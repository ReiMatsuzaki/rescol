#!/bin/sh
#PBS -l nodes=1:ppn=1,mem=1gb
#PBS -N L2Num40o4_R14.sh
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
time
python run.py
time
 
 
