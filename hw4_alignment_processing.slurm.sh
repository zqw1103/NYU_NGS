#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=44GB
#SBATCH --job-name=align_processing
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



#module load samtools/intel/1.9
# samtools view -bh $sample.sam > $sample.bam

module load picard/2.17.11

java -Xmx42g -XX:ParallelGCThreads=1 -jar "${PICARD_JAR}" SortSam \
INPUT=$sample.bam \
OUTPUT=$sample.sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR="${SLURM_JOBTMP}" \
MAX_RECORDS_IN_RAM=10000000 \
VALIDATION_STRINGENCY=LENIENT

echo _ESTATUS_ [ SortSam $sample ]: $?

module load samtools/intel/1.9

samtools index $sample.sorted.bam

java -Xmx42g -XX:ParallelGCThreads=1 -jar "${PICARD_JAR}" MarkDuplicates \
INPUT=$sample.sorted.bam \
OUTPUT=$sample.sorted.markdups.bam \
METRICS_FILE=$sample.metrics.txt \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT \
READ_NAME_REGEX=null \
TMP_DIR="${SLURM_JOBTMP}" \
MAX_RECORDS_IN_RAM=10000000 \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

echo _ESTATUS_ [ MarkDuplicates $sample ]: $?

samtools index $sample.sorted.markdups.bam

rm $sample.sorted.bam
rm $sample.sorted.bam.bai

######


echo _END_ [ alignment_processing for $sample ]: $(date)
