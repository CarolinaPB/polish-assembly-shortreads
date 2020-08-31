configfile: "config.yaml"
import os

ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])
BAM = os.path.join(config["DATADIR"], config["BAM"])

regions_list = [line.rstrip() for line in open(config["REGIONS"])]
regions_list = regions_list[:-1]
regions_list_replace = [el.replace("|", "__") for el in regions_list]


rule all:
    input:
        expand("variant_calling/{region}.vcf.gz", region=regions_list),
        # "results/bcftools_stats.txt"
        # expand("variant_calling/sorted/{region}.sorted.vcf.gz", region=regions_list)
        expand("variant_calling/{region}.sorted.indexed.txt", region=regions_list)

rule longshot:
    input:
        reference= ASSEMBLY, 
        bam=os.path.join(config["DATADIR"], config["BAM"])
    output:
        "variant_calling/{region}.vcf.gz" #make temporary
    shell:
        'longshot -r "{wildcards.region}" --bam {input.bam} --ref {input.reference} --out "{output}" --min_alt_count 20 --min_allele_qual 10 --min_alt_frac 0.2'

rule sort_vcf:
    input:
        rules.longshot.output
    output:
        "variant_calling/{region}.sorted.vcf.gz"
    shell:
        "module load bcftools && bcftools sort '{input}' > '{output}'"

rule index_vcf:
    input: 
        rules.sort_vcf.output
    output:
        temp("variant_calling/{region}.sorted.indexed.txt")
    shell:
        "module load bcftools && tabix -p vcf {input} && echo '{wildcards.region} indexed' > {output}"

rule merge_vcf:
    input:
        is_sorted=rules.sort_vcf.output,
        indexed=rules.index_vcf.output
    output:
        "variant_calling/var.merged.vcf.gz"
    shell:
        "module load vcftools && vcf-concat {input.is_sorted} -Oz -o {output}"

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
