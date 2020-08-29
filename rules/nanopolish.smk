rule nanopolish_index:
    input:
        nanopore_reads     = lambda wildcards: get_nanopore_reads(subsampled_reads, wildcards.sample, wildcards.subsampling, wildcards.coverage),
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
        nanopolish index -d {input.fast5_dir} -s {input.sequencing_summary} {output.linked_reads} 2>{log}
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
    threads: 16
    log: "logs/run_nanopolish/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 40000 * attempt
    singularity:
        nanopolish_container
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.ref} {input.nanopore_reads} 2> {log} | \
            samtools sort -o reads.sorted.bam -T reads.tmp 2>{log}
        samtools index reads.sorted.bam 2>{log}
        
        mkdir -p nanopolish_out/vcf        
        nanopolish_makerange.py {input.ref} | \
            parallel --results nanopolish_out -P 8 \
            nanopolish variants \
              -t 2 \
              -w {{1}} \
              --reads {input.nanopore_reads} \
              --bam reads.sorted.bam \
              --genome {input.ref} \
              -o nanopolish_out/vcf/nanopolish.{{1}}.vcf \
              -q dam,dcm \
              --ploidy 1
        cat nanopolish_out/vcf/nanopolish.*.vcf > {output.vcf}

        cp {input.ref} {output.ref}
        """

if "run_nanopolish_locally" in config and config["run_nanopolish_locally"]==True:
    print("Running nanopolish locally...")
    localrules: nanopolish_index, run_nanopolish
else:
    print("Submitting nanopolish jobs...")
