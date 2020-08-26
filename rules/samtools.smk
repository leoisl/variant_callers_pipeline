rule copy_ref:
    input:
        original_ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        linked_ref = output_folder+"/samtools/illumina/{coverage}x/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.ref.fa",
    threads: 1
    shell:
        "cp {input.original_ref} {output.linked_ref}"
localrules: copy_ref

rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    singularity: samtools_container
    shell:
        "bwa index {input.fasta} >{log} 2>&1"


rule faidx:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.fai"
    threads: 1
    log: "{fasta}.faidx.log"
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    singularity: samtools_container
    shell:
        "samtools faidx {input.fasta} >{log} 2>&1"


rule bwa_mem_map_reads_to_ref:
    input:
        illumina_reads_1 = lambda wildcards: get_illumina_reads(subsampled_reads, wildcards.sample, 1, wildcards.subsampling, wildcards.coverage),
        illumina_reads_2 = lambda wildcards: get_illumina_reads(subsampled_reads, wildcards.sample, 2, wildcards.subsampling, wildcards.coverage),
        ref = rules.copy_ref.output.linked_ref,
        ref_index = rules.copy_ref.output.linked_ref+".amb",
    output:
        bam = output_folder+"/samtools/illumina/{coverage}x/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.bam",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    log:
        "logs/bwa_mem_map_reads_to_ref/samtools/illumina/{coverage}x/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.bam.log"
    singularity: samtools_container
    shell:
        "bwa mem -t {threads} {input.ref} {input.illumina_reads_1} {input.illumina_reads_2} 2>{log} | samtools sort -@{threads} -o {output.bam} - 2>{log}"


rule samtools_mpileup_bcftools_call:
    input:
        ref = rules.copy_ref.output.linked_ref,
        ref_index = rules.copy_ref.output.linked_ref+".fai",
        bam = rules.bwa_mem_map_reads_to_ref.output.bam,
    output:
        vcf = output_folder+"/samtools/illumina/{coverage}x/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.vcf",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    log:
        "logs/samtools_mpileup_bcftools_call/samtools/illumina/{coverage}x/{subsampling}/{sample}/samtools_{sample}_AND_{reference}.vcf.log"
    singularity: samtools_container
    shell:
        "samtools mpileup -ugf {input.ref} {input.bam} 2>{log} | bcftools call --ploidy 1 -v -m -f GQ,GP -O v -o {output.vcf} 2>{log}"
