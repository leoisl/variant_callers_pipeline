rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    singularity: bwa_container
    shell:
        "bwa index {input.fasta} > {log} 2>&1"


rule faidx:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.fai"
    threads: 1
    log: "{fasta}.faidx.log"
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    singularity: samtools_container
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"


def get_ref_file_given_id(reference_id):
    return references_df.query("reference_id == @reference_id")["uncompressed_file"].to_list()[0]

rule map_reads_to_ref:
    input:
        reads_1=sample_data_dir+"/{sample}/{sample}.{coverage}.{subsampling}.illumina.1.fastq",
        reads_2=sample_data_dir+"/{sample}/{sample}.{coverage}.{subsampling}.illumina.2.fastq",
        ref = lambda wildcards: get_ref_file_given_id(wildcards.reference),
        ref_index = lambda wildcards: f"{get_ref_file_given_id(wildcards.reference)}.amb",
    output:
        bam = output_folder+"/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.bam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/map_reads_to_ref/{technology}.{coverage}.{subsampling}.{sample}.{reference}.log"
    singularity: bwa_container
    shell:
        "bwa mem -t {threads} {input.ref} {input.reads_1} {input.reads_2} | samtools sort -@{threads} -o {output.bam} -"


rule samtools_mpileup:
    input:
        ref=lambda wildcards: get_ref_file_given_id(wildcards.reference),
        ref_index=lambda wildcards: f"{get_ref_file_given_id(wildcards.reference)}.fai",
        bam=rules.map_reads_to_ref.output.bam,
    output:
        pileup = output_folder+"/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.pileup"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/call_variants_with_samtools/{technology}.{coverage}.{subsampling}.{sample}.{reference}.log"
    shell:
        "samtools mpileup -ugf {input.ref} {input.bam} > {output.pileup}"

rule bcftools_call:
    input:
        ref=lambda wildcards: get_ref_file_given_id(wildcards.reference),
        pileup = rules.samtools_mpileup.output.pileup
    output:
        ref_for_pipeline = output_folder + "/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.ref.fa",
        vcf = output_folder + "/{technology}/{coverage}/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.vcf",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/call_variants_with_samtools/{technology}.{coverage}.{subsampling}.{sample}.{reference}.log"
    run:
        shell("ln -s {input.ref} {output.ref_for_pipeline}")
        shell("cat {input.pileup} | bcftools call --ploidy 1 -v -m -f GQ,GP -O v -o {output.vcf}")
