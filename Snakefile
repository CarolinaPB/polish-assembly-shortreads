configfile: "config.yaml"
import os

# workdir: config["OUTDIR"]

ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])
BAM = os.path.join(config["DATADIR"], config["BAM"])

regions_list = [line.rstrip() for line in open(config["REGIONS"])]
regions_list = regions_list[:-1] # remove last line - contains asterisk
regions_list_replace = [el.replace("|", "__") for el in regions_list]

nlines = 0
with open(config["REGIONS"], "r") as infile:
    for line in infile:
        nlines +=1

# nfiles = int(round(nlines/3))
nfiles = nlines

rule all:
    input:
        # "variant_calling_sorted/var.merged.indexed.sorted.vcf.gz",
        "results/bcftools_stats.txt",
        expand("pilon/region_{n}.fasta", n=["%.3d" % i for i in range(nfiles)])

rule longshot:
    input:
        reference= ASSEMBLY, 
        bam=os.path.join(config["DATADIR"], config["BAM"])
    output:
        "variant_calling/{region}.vcf" # make temporary
    message:
        "Rule {rule} processing"
    shell:
        'longshot -r "{wildcards.region}" --bam {input.bam} --ref {input.reference} --out "{output}" --min_alt_count 20 --min_allele_qual 10 --min_alt_frac 0.2'

rule zip_vcf:
    input:
        rules.longshot.output
    output:
        temp("variant_calling/{region}.vcf.gz")
    message:
        "Rule {rule} processing"
    shell:
        "module load samtools && bgzip -c '{input}'> '{output}'"

rule index_vcf:
    input: 
        rules.zip_vcf.output
    output:
        temp("variant_calling/{region}.vcf.gz.tbi") # make temporary
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && tabix -p vcf '{input}'"

rule sort_vcf:
    input:
        files = rules.zip_vcf.output,
        index = rules.index_vcf.output
    output:
        temp("variant_calling_sorted/{region}.indexed.sorted.vcf.gz")
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && bcftools sort -Oz '{input.files}' > '{output}'"

rule merge_vcf:
    input:
        # rules.sort_vcf.output make temporary
        expand("variant_calling_sorted/{region}.indexed.sorted.vcf.gz", region = regions_list)
    output:
        temp("variant_calling_sorted/merged/var.merged.vcf.gz")
    message:
        "Rule {rule} processing"
    shell:
        "module load vcftools samtools && vcf-concat {input:q} | bgzip -c > {output}"

rule index_merged_vcf:
    input: 
        rules.merge_vcf.output
    output:
        "variant_calling_sorted/merged/var.merged.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && tabix -p vcf {input}"
        
rule sort_merged_vcf:
    input:
        idx = rules.index_merged_vcf.output,
        res = rules.merge_vcf.output
    output:
        "variant_calling_sorted/merged/var.merged.sorted.vcf.gz"
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && bcftools sort -Oz {input.res} > {output}"

rule index_sorted_vcf:
    input: 
        rules.sort_merged_vcf.output
    output:
        "variant_calling_sorted/merged/var.merged.sorted.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && tabix -p vcf {input}"


rule bcftools_stats:
    input:
        sort = rules.sort_merged_vcf.output,
        index = rules.index_sorted_vcf.output
    output:
        "results/bcftools_stats.txt"
    message:
        "Rule {rule} processing"
    shell:
        "module load bcftools && bcftools stats {input.sort} > {output}"



rule split_regions:
    input:
        config["REGIONS"]
    output:
        # dynamic("regions/region_{n}")
        # directory("regions/"),
        regions = "regions/region_{n}" # creates many unnecessary jobs. It's supposed to be done in one job. If make it so that it's only one job then the next rule will not recognize the output
    params:
        prefix = "region_"
    message:
        "Rule {rule} processing"
    shell:
        # "split -l 3 -d {input} regions/{params.prefix}"
        "split -l 1 -d -a 3 {input} regions/{params.prefix}"


rule pilon:
    input:
        reference= ASSEMBLY, 
        # bam=os.path.join(config["DATADIR"], config["BAM"]),
        short_bam = config["SHORT_READS_BAM"],
        regions = "regions/region_{n}",
        
    output:
        "pilon/region_{n}.fasta"
    params:
        outdir = "pilon/",
        depth=0.8,
    message:
        "Rule {rule} processing"
    shell:
        """
        pilon -Xmx150G --genome {input.reference} --bam {input.short_bam} --outdir {params.outdir} --output region_{wildcards.n} --vcf --diploid --mindepth {params.depth} --targets {input.regions} --threads 30 --fix "bases" 
        """