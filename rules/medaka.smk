rule run_medaka:
    input:
        nanopore_reads = lambda wildcards: get_nanopore_reads(samples, wildcards.sample, wildcards.subsampling, wildcards.coverage),
        ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        vcf = output_folder+"/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.vcf",
        ref = output_folder+"/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.ref.fa",
    shadow: "shallow"
    threads: 4
    log: "logs/run_medaka/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 20000 * attempt
    singularity:
        medaka_container
    shell:
        """
        # Following ARTIC pipeline from https://github.com/nanoporetech/medaka/issues/173#issuecomment-642566213 :
        minimap2 -t {threads} -ax map-ont {input.ref} {input.nanopore_reads} 2>{log} | samtools view -b - | samtools sort -@ {threads} -o nanopore_mapped.sorted.bam 2>{log}
        samtools index nanopore_mapped.sorted.bam 2>{log}
        medaka consensus --model r941_min_high_g344 nanopore_mapped.sorted.bam medaka_consensus_out.hdf 2>{log}
        medaka variant --verbose {input.ref} medaka_consensus_out.hdf medaka_variants.vcf 2>{log}
        
        # this is commented out because we have a bug (see https://github.com/iqbal-lab/pandora1_paper/issues/121#issuecomment-661879246 )        
        # medaka tools annotate medaka_variants.vcf {input.ref} nanopore_mapped.sorted.bam medaka_variants.annotated.vcf 2>{log}
        
        cp medaka_variants.vcf {output.vcf}
        cp {input.ref} {output.ref}
        """
