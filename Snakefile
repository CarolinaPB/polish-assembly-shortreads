configfile: "config.yaml"
import os

workdir: config["OUTDIR"]

#ASSEMBLY=os.path.join(config["DATADIR"], config["assembly"])
ASSEMBLY=config["assembly"]
BAM = os.path.join(config["DATADIR"], config["BAM"])

# regions_list = [line.rstrip() for line in open(config["REGIONS"])]
# regions_list = regions_list[:-1] # remove last line - contains asterisk
# regions_list_replace = [el.replace("|", "__") for el in regions_list]

nlines = 0
with open(config["REGIONS"], "r") as infile:
    for line in infile:
        nlines +=1

# nfiles = int(round(nlines/3))
nfiles = nlines

rule all:
    input:
        # "variant_calling_sorted/var.merged.indexed.sorted.vcf.gz",
        # "results/bcftools_stats.txt",
        # expand("pilon/region_{n}.{ext}", n=["%.3d" % i for i in range(nfiles)], ext=["fasta", "changes"]),
        "results/pilon.fasta",
        "results/pilon.changes",
        "results/pilon.amb"

rule split_regions:
    output:
        regions = expand(temp("regions/region_{n}"), n=["%.3d" % i for i in range(nfiles)])
    params:
        prefix = "region_",
        regions = config["REGIONS"]
    message:
        "Rule {rule} processing"
    shell:
        "split -l 1 -d -a 3 {params.regions} regions/{params.prefix}"


rule pilon:
    input:
        reference= ASSEMBLY, 
        short_bam = config["SHORT_READS_BAM"],
        regions = "regions/region_{n}"
        
    output: 
        temp(expand("pilon/region_{{n}}.{ext}", ext=["fasta", "changes"])) #make temporary
    params:
        outdir = "pilon/",
        depth=0.7,
    message:
        "Rule {rule} processing"
    log:
        "logs_rules/pilon/region_{n}.log"
    shell:
        """
        pilon -Xmx150G --genome {input.reference} --bam {input.short_bam} --outdir {params.outdir} --output region_{wildcards.n} --diploid --mindepth {params.depth} --targets {input.regions} --threads 30 --fix bases --changes > {log} 
        """

rule concat_pilon_fasta_changes:
    input:
        expand("pilon/region_{n}.{{ext}}", n=["%.3d" % i for i in range(nfiles)]),
        # expand("pilon/region_{n}.changes", n=["%.3d" % i for i in range(nfiles)]),
    output:
        "results/pilon.{ext}"
    message:
        "rule {rule} processing"
    shell:
        "cat {input} > {output}"

rule bwa_index:
    input: 
        # rules.concat_pilon_fasta.output
        "results/pilon.fasta"
    output:
        "results/pilon.amb"
    shell:
        "module load bwa && bwa index {input}"