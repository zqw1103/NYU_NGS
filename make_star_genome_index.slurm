#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --mem=115GB
#SBATCH --job-name=genomeGenerate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jmf11@nyu.edu

module purge

module load star/intel/2.7.6a

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir tmp_star \
--genomeFastaFiles Pdac_Barhee_chr_unan_cp_180126.fa \
--genomeSAindexNbases 13 \
--sjdbGTFfile Pdac_Barhee_chr_unan_180126_maker_HC.STAR_htseq.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene gene_id \
--sjdbOverhang 75

echo _ESTATUS_ [ genomeGenerate ]: $?
echo _END_ [ make_star_genome_index.slurm ]: $(date)
