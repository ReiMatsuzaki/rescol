#!/bin/csh -f
#PBS -l select=ncpus=1:mpiprocs=1:ompthreads=1:jobtype=large
#PBS -l walltime=72:00:00
cd ${PBS_O_WORKDIR}
python run.py

