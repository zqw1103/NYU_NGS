#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=bwamem_array
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<your email>
#SBATCH --array=1-30

module purge

# Please replace <your email> in the --mail-user directive above

# Define ref variable with path to normalized hg18 reference fasta
ref=

# Path to 3-column (tab-delimited) table with sample name, fastq 1 file name, and fastq 2 file name

table=/scratch/work/courses/BI7653/hw3.2021/fastqs.processed/hw3_fastqs.processed.txt
fqdir=/scratch/work/courses/BI7653/hw3.2021/fastqs.processed

# The following code defines sample, fq1 and fq2 variables for current array index
# note: SLURM_ARRAY_TASK_ID environmental variable will contain a single value corresponding to the current array index

line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample="$(printf "%s" "${line}" | cut -f1)"
fq1="$(printf "%s" "${line}" | cut -f2)"
fq2="$(printf "%s" "${line}" | cut -f3)"

# Print to standard out the array index and the sample name

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

# Make a directory for the sample and cd to it

mkdir $sample
cd $sample

# Load the bwa mem module

module load bwa/intel/0.7.17

bwa mem \
-M \
-t $SLURM_CPUS_PER_TASK \
-R "@RG\tID:${sample}.id\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}.lb" \
"${ref}" \
$fqdir/$fq1 \
$fqdir/$fq2 > $sample.sam

echo _ESTATUS_ [ bwa mem for $sample ]: $?

echo _END_ [ hw3_bwamem.slurm for $sample ]: $(date)
