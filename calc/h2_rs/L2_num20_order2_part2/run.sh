#!/bin/sh
#PBS -l nodes=1:ppn=1,mem=1gb
#PBS -N L2Num20o2.sh
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
python run.py
 
 
