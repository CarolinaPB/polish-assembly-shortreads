#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --exclude=fat001,fat002,fat101,fat100
#SBATCH --mem=16000


REFERENCE=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/HiC_assembly/Mgal6.fa
BAM=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/realign_pacbio/pb_mapped_Mgal6.sorted.bam
OUTPUT=/home/WUR/moiti001/polish-assembly-shortreads/variant_calling/bash_try_parallel2/var.vcf.gz

module load samtools vcflib/gcc/64/0.00.2019.07.10

conda init bash
conda activate freebayestest

time freebayes-parallel <(fasta_generate_regions.py ${REFERENCE}.fai 100000) 30 \
-f ${REFERENCE} --use-best-n-alleles 4 --min-base-quality 10 \
--min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 \
--min-alternate-count 2 --bam $BAM | vcffilter -f 'QUAL > 20'\
 | bgzip -c > $OUTPUT
