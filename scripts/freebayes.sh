#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --exclude=fat001,fat002,fat101,fat100
#SBATCH --mem=16000


REFERENCE=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/HiC_assembly/Mgal6.fa
BAM=/lustre/backup/WUR/ABGC/shared/ABGC_Projects/Turkey_Assembly/realign_pacbio/pb_mapped_Mgal6.sorted.bam
OUTPUT=/home/WUR/moiti001/polish-assembly-shortreads/variant_calling/var.vcf.gz

module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10

time freebayes -f $REFERENCE -r "Mgal6_Chr7|Scaffold_8__1_contigs__length_36428186" --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam $BAM | vcffilter -f 'QUAL > 20' | bgzip -c > $OUTPUT