rule run_nanopolish:
    input:
        nanopore_reads     = lambda wildcards: get_nanopore_reads(subsampled_reads, wildcards.sample, wildcards.subsampling, wildcards.coverage),
        # Obs: fine to give a fastq and a superset of fast5s: https://github.com/jts/nanopolish/issues/443
        fast5_dir          = lambda wildcards: get_fast5_dir(fast5s_df, wildcards.sample),
        sequencing_summary = lambda wildcards: get_sequencing_summary(fast5s_df, wildcards.sample),
        ref = lambda wildcards: get_uncompressed_reference(references, wildcards.reference),
    output:
        copied_reads = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq",
        index = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index",
        fai = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.fai",
        gzi = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.gzi",
        readdb = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.reads.fastq.index.readdb",
        vcf = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.vcf",
        ref = output_folder+"/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.ref.fa",
    shadow: "shallow"
    threads: 16
    params:
        nb_of_processes = 8
    log: "logs/run_nanopolish/nanopolish/nanopore/{coverage}x/{subsampling}/{sample}/nanopolish_{sample}_AND_{reference}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt
    singularity:
        nanopolish_container
    shell:
        """
        cp {input.nanopore_reads} {output.copied_reads}
        nanopolish index -d {input.fast5_dir} -s {input.sequencing_summary} {output.copied_reads} 2>{log}

        minimap2 -ax map-ont -t {threads} {input.ref} {output.copied_reads} 2> {log} | \
            samtools sort -o reads.sorted.bam -T reads.tmp 2>{log}
        samtools index reads.sorted.bam 2>{log}
        
        mkdir -p nanopolish_out/vcf        
        nanopolish_makerange.py {input.ref} | \
            parallel --results nanopolish_out -P {params.nb_of_processes} \
            nanopolish variants \
              -t $(({threads} / {params.nb_of_processes})) \
              -w {{1}} \
              --reads {output.copied_reads} \
              --bam reads.sorted.bam \
              --genome {input.ref} \
              -o nanopolish_out/vcf/nanopolish.{{1}}.vcf \
              -q dam,dcm \
              --ploidy 1 \
              2>{log}
        
        # concat vcfs
        cat `ls -1 nanopolish_out/vcf/nanopolish.*.vcf | head -n 1` | grep "^#" > nanopolish_out/vcf_header
        cat nanopolish_out/vcf/nanopolish.*.vcf | grep -v "^#" > nanopolish_out/vcf_content
        cat nanopolish_out/vcf_header nanopolish_out/vcf_content > {output.vcf} 
        
        cp {input.ref} {output.ref}
        """

if "run_nanopolish_locally" in config and config["run_nanopolish_locally"]==True:
    print("Running nanopolish locally...")
    localrules: nanopolish_index, run_nanopolish
else:
    print("Submitting nanopolish jobs...")
