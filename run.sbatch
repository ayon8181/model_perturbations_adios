#!/bin/bash
#
#SBATCH --job-name=1D_ref
#SBATCH --nodes=8
#SBATCH --ntasks=384
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
## project
#SBATCH -A EAR21003
## job queue: normal, development, rtx, rtx-dev
#SBATCH -p development
#SBATCH -o 1D_ref.out

##################################################

NPROC=384

##################################################

echo "Module list"
module list

echo "On which nodes it executes" 
echo $SLURM_JOB_NODELIST 

# begin mesher
echo
echo "running mesher: `date`"
echo "directory: `pwd`"
echo

ibrun -n $NPROC ./bin/xmeshfem3D
#/home1/09038/ayon8181/model_perturbation_adios/make all
#cp /home1/09038/ayon8181/model_perturbation_adios/bin/* ./bin/
#echo "Creating perturbation model"
#ibrun -n $NPROC ./bin/xmodel_pert 10 10 2000 1000 /scratch1/09038/ayon8181/s40rts_deep/DATABASES_MPI /scratch1/09038/ayon8181/s40rts_deep/DATABASES_MPI '1 1 1 1 1 1 1'
#echo "Creating Elastic Perturbations"
#ibrun -n $NPROC ./bin/xmodel_add DATABASES_MPI DATABASES_MPI DATABASES_MPI 0.01
#echo "Creating Qmu perturbation"
#ibrun -n $NPROC ./bin/xmodel_add_qmu DATABASES_MPI DATABASES_MPI DATABASES_MPI 10

## checks exit code
#if [[ $? -ne 0 ]]; then exit 1; fi

echo "done: `date`"
# end mesher

# begin solver
#echo
echo "running solver: `date`"
echo
 
#ibrun -n $NPROC ./bin/xspecfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done: `date`"
echo
## end solver
