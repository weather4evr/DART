#!/bin/csh

#
#BSUB -n 81
##BSUB -n 16
#BSUB -J convert
#BSUB -o output.driver
#BSUB -e output.driver
#BSUB -q regular
#BSUB -P P64000510
#BSUB -W 15
#BSUB -R "span[ptile=16]"

#PBS -S /bin/csh
#PBS -N test
#PBS -A NMMM0021
#PBS -l walltime=10:00
#PBS -q economy
#PBS -o ./output_file
#PBS -j oe 
#PBS -k eod 
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -m n    
#PBS -M schwartz@ucar.edu
#PBS -V 
#

cd $PWD

rm -rf core*
rm -f obs_seq.00??
mpiexec_mpt ./jedi_to_dart
