configfile: "config.yaml"
import os

workdir: config["OUTDIR"]

include: "/lustre/nobackup/WUR/ABGC/moiti001/snakemake-rules/create_file_log.smk"

pipeline = "polish-assembly-shortreads"


ASSEMBLY=config["assembly"]
BAM = config["BAM"]

localrules: split_regions, get_corrected_summary, create_regions, create_file_log

rule all:
    input:
        files_log,
        expand("results/pilon.{ext}", ext=["fasta", "changes"]),
        "results/pilon.fasta.amb",
        "results/summary_corrected_sites.txt",
        files_log

rule create_regions:
    input:
        BAM
    output:
        "chr_names.txt"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load samtools && samtools idxstats {input} | cut -f 1 > {output}
        awk '!/*/' {output} > temp && mv temp {output}
        """

rule split_regions:
    input:
        "chr_names.txt"
    output:
        regions = temp(expand('regions/region_{n}', n=["%.3d" % i for i in range(config["nregions"])]))
    params:
        prefix = "region_"
        # regions = config["REGIONS"]
    message:
        "Rule {rule} processing"
    shell:
        "split -l 1 -d -a 3 {input} regions/{params.prefix}"

rule pilon:
    input:
        reference= ASSEMBLY,
        bam = BAM,
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
        pilon -Xmx150G --genome {input.reference} --bam {input.bam} --outdir {params.outdir} --output region_{wildcards.n} --diploid --mindepth {params.depth} --targets {input.regions} --threads 30 --fix bases --changes > {log} 
        """

rule concat_pilon_fasta_changes:
    input:
        expand("pilon/region_{n}.{{ext}}", n=["%.3d" % i for i in range(config["nregions"])]),
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
        "results/pilon.fasta.amb"
    shell:
        "module load bwa && bwa index {input}"

rule get_corrected_summary:
    input:
        "results/pilon.fasta"
    params:
        logs = "logs_rules/pilon"
    output:
        "results/summary_corrected_sites.txt"
    message:
        "Rule {rule} processing"
    run:
        snps = 0
        insertions=0
        deletions=0
        log_dir=params.logs

        for file in os.listdir(log_dir):
            with open(os.path.join(log_dir,file), "r") as infile:
                for line in infile:
                    if line.startswith("Corrected"):
                        extracted = [int(s) for s in line.split() if s.isdigit()]
                        snps += extracted[0]
                        insertions += extracted[1]
                        deletions += extracted[3]

        with open(output[0], "w") as out:
            out.write("Corrected\n")
            out.write(f"snps: {snps}\n")
            out.write(f"insertions: {insertions}\n")
            out.write(f"deletions: {deletions}\n")
