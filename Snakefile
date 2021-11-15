################################################
#ATACSeq snakemake WF
# Author: DIOP Khadidiatou
# HPC PROJECT
################################################




configfile: "config/config.yaml",

rule all:
    input:
        expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
        expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
        expand("results/trimming/{echan}_1.clean.fastq", echan = config["echan"]),
        expand("results/trimming/{echan}_2.clean.fastq", echan = config["echan"]),
        expand("results/fastqc_post_trim/{sample}.clean_fastqc.zip", sample = config["samples"]),
        expand("results/fastqc_post_trim/{sample}.clean_fastqc.html", sample = config["samples"]),
        expand("results/bowtie2/{echan}_mapped.bam", echan = config["echan"]),
        expand("results/picard/{echan}_mapped_clean.bam", echan = config["echan"]),
        expand("results/picard/{echan}_mapped_net.txt", echan = config["echan"]),
        expand("results/deeptools/PlotCoverage.pdf"),
        expand("results/deeptools/resultsMultiBamSummary.npz"),
        expand("results/deeptools/heatmap_SpearmanCorr_readCounts.png"),
        expand("results/macs/{echan}_mapped_clean_macs_peaks.xls", echan = config["echan"]),
        expand("results/macs/{echan}_mapped_clean_macs_summits.bed", echan = config["echan"]),
        expand("results/macs/{echan}_mapped_clean_macs_peaks.narrowPeak", echan = config["echan"]),
        expand("results/macs/{echan}_mapped_clean_macs_model.r", echan = config["echan"]),
        expand("results/bedtools/ss_50k_{id}_0h_24h_common_peaks.bed", id=["1","2","3"]),
        expand("results/bedtools/ss_50k_{id}_0h_unique_peaks.bed", id=["1","2","3"]),
        expand("results/bedtools/ss_50k_{id}_24h_unique_peaks.bed", id=["1","2","3"])

rule unzip:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("tmp/{sample}.fastq")
    shell:
        """
        mkdir -p tmp
        gunzip -c {input} > {output}
        """

rule fastqc_init:
    input:
        "tmp/{sample}.fastq"
    output:
        html="results/fastqc_init/{sample}_fastqc.html",
        zip="results/fastqc_init/{sample}_fastqc.zip"
    conda:
       "env.yaml"
    threads: 2
    shell:
        """
        mkdir -p results/fastqc_init
        fastqc {input} -o "results/fastqc_init" -t {threads}
        """


rule cutadapt:
    input:
        read="tmp/{echan}_1.fastq",
        read2="tmp/{echan}_2.fastq"
    output:
        clean_R1="results/trimming/{echan}_1.clean.fastq",
        clean_R2="results/trimming/{echan}_2.clean.fastq"
    threads: 8
    conda:
       "env.yaml"

    shell:
        """
        mkdir -p results/trimming
        cutadapt -j {threads} -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o {output.clean_R1} -p {output.clean_R2} {input.read} {input.read2}
        """

rule fastqc_post_trim:
    input:
        "results/trimming/{sample}.clean.fastq"
    output:
        html="results/fastqc_post_trim/{sample}.clean_fastqc.html",
        zip="results/fastqc_post_trim/{sample}.clean_fastqc.zip"
    conda:
        "env.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/fastqc_post_trim
        fastqc {input} -o "results/fastqc_post_trim" -t {threads}
        """


rule bowtie2:
    input:
        R1="results/trimming/{echan}_1.clean.fastq",
        R2="results/trimming/{echan}_2.clean.fastq"
    output:
        bam="results/bowtie2/{echan}_mapped.bam"
    conda:
         "env.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/bowtie2
        bowtie2  --very-sensitive -p {threads} -x data/mydatalocal/atacseq/bowtie2/all -1 {input.R1} -2 {input.R2} |  samtools view -q 2 -bS  -  |  samtools sort - -o {output.bam}
        """

rule picard:
    input:
        map="results/bowtie2/{echan}_mapped.bam"
    output:
        bamnet="results/picard/{echan}_mapped_clean.bam",
        net_txt="results/picard/{echan}_mapped_net.txt"

    conda:
        "env.yaml"
    threads: 6
    shell:
        """
        mkdir -p results/picard
        picard MarkDuplicates -I {input.map} -O {output.bamnet} -M {output.net_txt} -REMOVE_DUPLICATES true
        samtools index -b {output.bamnet}
        """

rule deeptools:
    params:
        bam="results/picard/*.bam"
    output:
        plot_cov="results/deeptools/PlotCoverage.pdf",
        summary="results/deeptools/resultsMultiBamSummary.npz",
        plot_corr="results/deeptools/heatmap_SpearmanCorr_readCounts.png"
    conda:
        "env.yaml"
    threads: 6
    shell:
        """
        mkdir -p results/deeptools
        plotCoverage --bamfiles {params.bam} --plotFile {output.plot_cov} --smartLabels --plotFileFormat pdf -p 6
        multiBamSummary bins --bamfiles {params.bam} -o {output.summary} -p 6
        plotCorrelation -in {output.summary} --corMethod spearman --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.plot_corr}
        """


rule macs2:
    input:
        pic="results/picard/{echan}_mapped_clean.bam"
    output:
        narrowpeaks="results/macs/{echan}_mapped_clean_macs_peaks.narrowPeak",
        bed="results/macs/{echan}_mapped_clean_macs_summits.bed",
        peaks_xls="results/macs/{echan}_mapped_clean_macs_peaks.xls",
        r="results/macs/{echan}_mapped_clean_macs_model.r"
    params:
        nom="{echan}_mapped_clean_macs",
        rep="results/macs"
    conda:
        "env.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/macs
        macs2 callpeak -t {input.pic} -f BAM -n {params.nom} --outdir {params.rep}
        """

rule bedtools:
    input:
        zero="results/macs/ss_50k_0h_R{id}_mapped_clean_macs_summits.bed",
        V="results/macs/ss_50k_24h_R{id}_mapped_clean_macs_summits.bed"
    output:
        commun="results/bedtools/ss_50k_{id}_0h_24h_common_peaks.bed",
        uniquezero="results/bedtools/ss_50k_{id}_0h_unique_peaks.bed",
        onlyV="results/bedtools/ss_50k_{id}_24h_unique_peaks.bed"
    conda:
        "env.yaml"
    threads: 6
    shell:
        """
        mkdir -p results/bedtools
        bedtools intersect -a {input.zero} -b {input.V} > {output.commun}
        bedtools intersect -v -a {input.zero} -b {input.V} > {output.uniquezero}
        bedtools intersect -v -a {input.V} -b {input.zero} > {output.onlyV}
        """

