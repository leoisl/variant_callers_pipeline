rule run_medaka:
    input:
        nanopore_reads = lambda wildcards: get_nanopore_reads(subsampled_reads, wildcards.sample, wildcards.subsampling, wildcards.coverage),
        ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        vcf = output_folder+"/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.vcf",
        ref = output_folder+"/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.ref.fa",
    shadow: "shallow"
    threads: 4
    log: "logs/run_medaka/medaka/nanopore/{coverage}x/{subsampling}/{sample}/medaka_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    singularity:
        medaka_container
    shell:
        """
        # first check if it is the pair of sample and ref that always fails due to a software bug
        if [ "{wildcards.sample}" == "Escherichia_coli_MSB1_7A" ] && [ "{wildcards.reference}" == "NZ_CP011134.1" ]
            then
            # medaka will bug on this input, create an empty VCF
            # Note: we still don't want to automate this because medaka can fail for other reasons that are not software bug,
            # but machine related or user error
            cat > {output.vcf} <<- EndOfMessage
##fileformat=VCFv4.1
##medaka_version=1.0.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
EndOfMessage
        else
            # run medaka normally
            # Following ARTIC pipeline from https://github.com/nanoporetech/medaka/issues/173#issuecomment-642566213 :
            minimap2 -t {threads} -ax map-ont {input.ref} {input.nanopore_reads} 2>{log} | samtools view -b - | samtools sort -@ {threads} -o nanopore_mapped.sorted.bam 2>{log}
            samtools index nanopore_mapped.sorted.bam 2>{log}
            medaka consensus --model r941_min_high_g344 nanopore_mapped.sorted.bam medaka_consensus_out.hdf 2>{log}
            medaka variant --verbose {input.ref} medaka_consensus_out.hdf medaka_variants.vcf 2>{log}
            
            # this is commented out because we have a bug (see https://github.com/iqbal-lab/pandora1_paper/issues/121#issuecomment-661879246 )        
            # medaka tools annotate medaka_variants.vcf {input.ref} nanopore_mapped.sorted.bam medaka_variants.annotated.vcf 2>{log}
            
            cp medaka_variants.vcf {output.vcf}
        fi

        cp {input.ref} {output.ref}
        """
