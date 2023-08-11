#!/bin/bash
## Run Single MCMC
#SBATCH --job-name=rungrid
## Allocation Definition
## On mox and ikt, the account and partition options should be the same.
#SBATCH --account=astro 
#SBATCH --partition=astro
## Resources
## Nodes
#SBATCH --nodes=1   
## Walltime (3 hours). Do not specify a walltime substantially more than your job needs.
#SBATCH --time=0-10:00:00
## Memory per node. It is important to specify the memory since the default memory is very small.
## For mox, --mem may be more than 100G depending on the memory of your nodes.
## For ikt, --mem may be 58G or more depending on the memory of your nodes.
## See above section on "Specifying memory" for choices for --mem.
#SBATCH --mem=100G
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/astro/tagordon/hydrosphere/h2o_vs_temp
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tagordon@uw.edu
## export all your environment variables to the batch job session
#SBATCH --export=ALL

module load anaconda3_5.3
python3 RunGrid.py 0.75 0.87 100 /gscratch/astro/tagordon/hydrosphere/h2o_vs_temp/m_05_x_15 
