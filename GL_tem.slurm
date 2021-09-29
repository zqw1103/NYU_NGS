#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=GL_template
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge
echo script begin: $(date)

module load angsd/intel/0.933

angsd -bam bamfiles_for_GLs.txt -P ${SLURM_CPUS_PER_TASK} -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -doGlf 2 -out GLs.gz

echo _ESTATUS_ [ GL ]: $?
echo _END_ [ GL_tem.slurm ]: $(date)
