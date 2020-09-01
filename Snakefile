configfile: "config.yaml"
import os

# workdir: config["OUTDIR"]

ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])
BAM = os.path.join(config["DATADIR"], config["BAM"])

regions_list = [line.rstrip() for line in open(config["REGIONS"])]
regions_list = regions_list[:-1] # remove last line - contains asterisk
regions_list_replace = [el.replace("|", "__") for el in regions_list]


rule all:
    input:
        # expand("variant_calling/{region}.vcf.gz", region=regions_list),
        # "results/bcftools_stats.txt"
        # expand("variant_calling/sorted/{region}.sorted.vcf.gz", region=regions_list)
        # expand("variant_calling_sorted/{region}.indexed.sorted.vcf.gz", region=regions_list),
        "variant_calling_sorted/var.merged.vcf.gz"
rule longshot:
    input:
        reference= ASSEMBLY, 
        bam=os.path.join(config["DATADIR"], config["BAM"])
    output:
        temp("variant_calling/{region}.vcf")
    shell:
        'longshot -r "{wildcards.region}" --bam {input.bam} --ref {input.reference} --out "{output}" --min_alt_count 20 --min_allele_qual 10 --min_alt_frac 0.2'

rule zip_vcf:
    input:
        rules.longshot.output
    output:
        temp("variant_calling/{region}.vcf.gz")
    shell:
        "module load samtools && bgzip -c '{input}'> '{output}'"

rule index_vcf:
    input: 
        rules.zip_vcf.output
    output:
        "variant_calling/{region}.vcf.gz.tbi"
    shell:
        "module load bcftools && tabix -p vcf '{input}'"

rule sort_vcf:
    input:
        files = rules.zip_vcf.output,
        index = rules.index_vcf.output
    output:
        temp("variant_calling_sorted/{region}.indexed.sorted.vcf.gz")
    shell:
        "module load bcftools && bcftools sort '{input.files}' > '{output}'"

rule merge_vcf:
    input:
        # rules.sort_vcf.output
        expand("variant_calling_sorted/{region}.indexed.sorted.vcf.gz", region = regions_list)
    output:
        "variant_calling_sorted/var.merged.vcf.gz"
    shell:
        "module load vcftools samtools && vcf-concat {input:q} | bgzip -c > {output}"

rule bcftools_stats:
    input:
        rules.merge_vcf.output
    output:
        "results/bcftools_stats.txt"
    shell:
        "module load bcftools && bcftools stats {input} > {output}"


# rule pilon:
#     input:
#         reference= ASSEMBLY, 
#         bam=os.path.join(config["DATADIR"], config["BAM"])
#     output:
#         directory("pilon/")
#     message:
#         "Rule {rule} processing"
#     params:
#         depth=0.9
#     shell:
#         "pilon --genome {input.reference} --bam {input.bam} --outdir {output} --vcf --diploid --mindepth {params.depth} --threads 10 --fix all"
#         # "pilon --genome $REFERENCE --bam $BAM --outdir $OUTPUT --vcf --diploid --frags <bam from short reads> --mindepth {params.depth}"
