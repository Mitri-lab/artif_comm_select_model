#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1 
#SBATCH --mem 2G 
#SBATCH --time 03:44:00

module load gcc python
cd /scratch/pguridif/210822
source /work/FAC/FBM/DMF/smitri/evomicrocomm/pablo/bin/activate
seed=$1
method=$2
python 210817_model_strains_poisson.py $seed $method
