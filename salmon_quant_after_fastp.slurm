#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=salmon_quant
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
#SBATCH --array=1-6

module purge

module load salmon/1.4.0

echo The array index is: ${SLURM_ARRAY_TASK_ID}

table=/scratch/zw2946/final.proj/processed_fastq_file.txt

fqfilepath="$(head -n ${SLURM_ARRAY_TASK_ID} "${table}" | tail -n 1)"
fq="$(basename $fqfilepath)"
sample="$(basename $fq .sE.fastq.gz)"

salmon_index_dir=/scratch/zw2946/final.proj/Homo_sapiens.GRCh38.cdna.all.normalized.fa_index

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

mkdir "${sample}"
cd "${sample}"

salmon quant -i ${salmon_index_dir} -l A -r $fqfilepath  --validateMappings --gcBias \
--threads ${SLURM_CPUS_PER_TASK} -o ${sample}.transcripts_quant

echo _ESTATUS_ [ salmon quant $sample ]: $?
echo _END_ [ salmon_quant_after_fastp.slurm ]: $(date)
