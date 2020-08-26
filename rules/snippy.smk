rule run_pe_snippy:
    input:
        illumina_reads_1 = lambda wildcards: get_illumina_reads(subsampled_reads, wildcards.sample, 1, wildcards.subsampling, wildcards.coverage),
        illumina_reads_2 = lambda wildcards: get_illumina_reads(subsampled_reads, wildcards.sample, 2, wildcards.subsampling, wildcards.coverage),
        ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        vcf = output_folder+"/snippy/illumina/{coverage}x/{subsampling}/{sample}/snippy_{sample}_AND_{reference}.vcf",
        ref = output_folder+"/snippy/illumina/{coverage}x/{subsampling}/{sample}/snippy_{sample}_AND_{reference}.ref.fa",
    shadow: "shallow"
    threads: 4
    log: "logs/run_pe_snippy/illumina/{coverage}x/{subsampling}/{sample}/snippy_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    singularity:
        snippy_container
    shell:
        """
        snippy --cpus {threads} --outdir snippy_outdir --reference {input.ref} --pe1 {input.illumina_reads_1} --pe2 {input.illumina_reads_2} >{log} 2>&1
        cp snippy_outdir/snps.filt.vcf {output.vcf}
        cp {input.ref} {output.ref}
        """
