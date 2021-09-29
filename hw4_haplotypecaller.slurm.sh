#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=hc_Sep23
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<your email>
#SBATCH --array=1-30

module purge

# Please replace <your email> in the --mail-user directive above

# Define ref variable with path to normalized hg18 reference fasta
ref=/scratch/work/courses/BI7653/hw3.2021/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa

# Path to 3-column (tab-delimited) table with sample name, fastq 1 file name, and fastq 2 file name

table=/scratch/work/courses/BI7653/hw3.2021/fastqs.processed/hw3_fastqs.processed.txt


# Define sample, fq1 and fq2 variables for current array index
# note: SLURM_ARRAY_TASK_ID environmental variable will contain a single value corresponding to the current array index

line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample="$(printf "%s" "${line}" | cut -f1)"

# Print to standard out the array index and the sample name

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample



cd $sample

module load gatk/4.0.2.1

gatk --java-options "-Xmx6G" HaplotypeCaller \
-R $ref \
-I $sample.sorted.markdups.bam \
-O $sample.g.vcf \
-ERC GVCF

echo _ESTATUS_ [ HaplotypeCaller $sample ]: $?


echo _END_ [ hw4_haplotypecaller.slurm $sample ]: $(date)
