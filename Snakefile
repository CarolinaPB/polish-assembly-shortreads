configfile: "config.yaml"
import os

ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])

rule all:
    input:
        # "variant_calling/var.vcf.gz"
        "variant_calling/var.vcf"

rule freebayes:
    input: 
        reference= ASSEMBLY, 
        bam=os.path.join(config["DATADIR"], config["BAM"])
    conda:
        "envs/freebayes.yaml"
    output: 
        # "variant_calling/var.vcf.gz"
        "variant_calling/var.vcf"
    message:
        "Rule {rule} processing"
    shell:
        # "module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output}"
        # "module load python/2.7.15 samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes-parallel <(fasta_generate_regions.py {input.reference}.fai 100000) 16 -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20'| bgzip -c > {output}"
        "module load python/2.7.15 samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes-parallel <(fasta_generate_regions.py {input.reference}.fai 100000) 16 -f {input.reference} --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} > {output}"

# use python 2 since vcflib has some dependencies on it