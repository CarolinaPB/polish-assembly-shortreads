configfile: "config.yaml"
import os

workdir: config["OUTDIR"]

ASSEMBLY=config["assembly"]
BAM = config["BAM"]


nlines = 0
with open(config["REGIONS"], "r") as infile:
    for line in infile:
        nlines +=1

# nfiles = int(round(nlines/3))
nfiles = nlines

localrules: all, split_regions, prepare_assembly

rule all:
    input:
        expand("results/pilon.{ext}", ext=["fasta", "changes", "amb"])

rule prepare_assembly:
# If it's the second round of polishing, the .fasta file will have |pilon at the end of the chr names.
# This will remove it so the rest of the analysis can work
    input:
        ASSEMBLY
    output:
        temp("data/assembly.fasta")
    message:
        "Rule {rule} processing"
    params:
        pattern = "|pilon"
    shell:
        """
        cp {input} {output}
        if grep -Fq {params.pattern:q} {output}
            then
                sed -i -e 's/|pilon//g' {output}
        fi
        """

rule split_regions:
    output:
        regions = temp(expand(temp("regions/region_{n}"), n=["%.3d" % i for i in range(nfiles)]))
    params:
        prefix = "region_",
        regions = config["REGIONS"]
    message:
        "Rule {rule} processing"
    shell:
        "split -l 1 -d -a 3 {params.regions} regions/{params.prefix}"

rule pilon:
    input:
        reference= rules.prepare_assembly.output,
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
    output:
        "results/pilon.{ext}"
    message:
        "rule {rule} processing"
    shell:
        "cat {input} > {output}"

rule bwa_index:
    input: 
        "results/pilon.fasta"
    output:
        "results/pilon.amb"
    shell:
        "module load bwa && bwa index {input}"

