configfile: "config.yaml"
import os

ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])

rule all:
    input:
        "variant_calling/var.vcf.gz"

rule freebayes:
    input: 
        reference= ASSEMBLY, 
        bam=os.path.join(config["DATADIR"], config["BAM"])
    output: 
        "variant_calling/var.vcf.gz"
    message:
        "Rule {rule} processing"
    shell:
        "module load freebayes samtools vcflib/gcc/64/0.00.2019.07.10 && freebayes -f {input.reference} -r 'Mgal6_Chr30|Scaffold_28__1_contigs__length_5171530' --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 --bam {input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output}"
