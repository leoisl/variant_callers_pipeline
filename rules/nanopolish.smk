rule nanopolish_index:
    input:
        nanopore_reads     = lambda wildcards: get_nanopore_reads(samples, wildcards.sample, wildcards.subsampling, wildcards.coverage),
        # Obs: fine to give a fastq and a superset of fast5s: https://github.com/jts/nanopolish/issues/443
        fast5_dir          = lambda wildcards: get_fast5_dir(fast5s_df, wildcards.sample),
        sequencing_summary = lambda wildcards: get_sequencing_summary(fast5s_df, wildcards.sample),
    output:
        linked_reads = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq",
        index = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index",
        fai = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.fai",
        gzi = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.gzi",
        readdb = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.readdb",
    threads: 1
    log: "logs/nanopolish_index/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity:
        nanopolish_container
    shell:
        """
        ln -s {input.nanopore_reads} {output.linked_reads}
        nanopolish index -d {input.fast5_dir} -s {input.sequencing_summary} {output.linked_reads}
        """


rule run_nanopolish:
    input:
        nanopore_reads = rules.nanopolish_index.output.linked_reads,
        nanopolish_index = rules.nanopolish_index.output.index,
        ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        vcf = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.vcf",
        ref = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.ref.fa",
    shadow: "shallow"
    threads: 32
    log: "logs/run_nanopolish/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 40000 * attempt
    singularity:
        nanopolish_container
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.ref} {input.nanopore_reads} | \
            samtools sort -o reads.sorted.bam -T reads.tmp
        samtools index reads.sorted.bam
        
        nanopolish variants \
          -t {threads} \
          --reads {input.nanopore_reads} \
          --bam reads.sorted.bam \
          --genome {input.ref} \
          -o {output.vcf} \
          -q dam,dcm \
          --ploidy 1
      
        cp {input.ref} {output.ref}
        """
