configfile: "config.yaml"

rule all:
        input:
            expand('lambda.{sample}_1_bismark_bt2_pe.bam', sample=config["samples"]),
            expand('{sample}.CpG_report.merged_CpG_evidence.cov', sample=config["samples"])
      
rule alignment:
    input:
        read1="{sample}_1.fq.gz",
        read2="{sample}_2.fq.gz",
        genome_dir="genome_folder/Bisulfite_Genome"
    output:
        "{sample}_1_bismark_bt2_pe.bam"
    shell:
        """
        /bin/Bismark-0.22.3/bismark --multicore 3 --genome genome_folder -1 {input.read1} -2 {input.read2}
        """

rule lambda_alignment:
    input:
        read1="{sample}_1.fq.gz",
        read2="{sample}_2.fq.gz",
        genome_dir="lambda_genome/Bisulfite_Genome"
    output:
        "lambda.{sample}_1_bismark_bt2_pe.bam"
    shell:
        """
        /bin/Bismark-0.22.3/bismark --multicore 3 --genome lambda_genome --prefix lambda -1 {input.read1} -2 {input.read2}
        """

rule deduplication:
    input:
        "{sample}_1_bismark_bt2_pe.bam"
    output:
        "{sample}_1_bismark_bt2_pe.deduplicated.bam"
    shell:
        """
        /bin/Bismark-0.22.3/deduplicate_bismark {input}
        """

rule meth_extraction:
    input:
        "{sample}_1_bismark_bt2_pe.deduplicated.bam"
    output:
        "{sample}_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    shell:
        """
        /bin/Bismark-0.22.3/bismark_methylation_extractor \
        --multicore 3 --scaffolds \
        --no_overlap --comprehensive --merge_non_CpG --bedgraph --report --cytosine_report \
        --genome_folder genome_folder {input}
        """

rule merge_cpgs:
    input:
        "{sample}_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    output:
        "{sample}.CpG_report.merged_CpG_evidence.cov"
    params:
        "{sample}"
    shell:
        """
        /bin/Bismark-0.22.3/coverage2cytosine \
        -o {params} --merge_CpGs \
        --genome_folder genome_folder {input}
        """
