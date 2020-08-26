import pandas as pd
from pathlib import Path

# ======================================================
# Utility functions
# ======================================================
def update_to_absolute_path_core(path_series):
    return path_series.apply(lambda path: str(Path(path).absolute()))
def update_to_absolute_path(df, columns):
    for column in columns:
        df[column] = update_to_absolute_path_core(df[column])
    return df

def get_illumina_reads(subsampled_reads, sample_name, first_or_second, subsampling, coverage):
    assert first_or_second in [1, 2]
    assert sample_name in subsampled_reads.sample_id.to_list()
    sample_path = subsampled_reads[subsampled_reads.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.{coverage}x.{subsampling}.illumina.{first_or_second}.fastq"

def get_nanopore_reads(subsampled_reads, sample_name, subsampling, coverage):
    assert sample_name in subsampled_reads.sample_id.to_list()
    sample_path = subsampled_reads[subsampled_reads.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.{coverage}x.{subsampling}.nanopore.fastq"

def get_uncompressed_reference(references, reference_id):
    assert reference_id in references.reference_id.to_list()
    uncompressed_reference_path = references[references.reference_id == reference_id]["uncompressed_file"].tolist()[0]
    return uncompressed_reference_path

def get_fast5_dir(fast5s_df, sample_id):
    assert sample_id in fast5s_df.sample_id.to_list()
    fast5_dir = fast5s_df[fast5s_df.sample_id == sample_id]["fast5_dir"].tolist()[0]
    return fast5_dir

def get_sequencing_summary(fast5s_df, sample_id):
    assert sample_id in fast5s_df.sample_id.to_list()
    sequencing_summary = fast5s_df[fast5s_df.sample_id == sample_id]["sequencing_summary"].tolist()[0]
    return sequencing_summary



# ======================================================
# Config files and vars
# ======================================================
configfile: "config.yaml"
output_folder = config['output_folder']
illumina_tools = config["illumina_tools"]
nanopore_tools = config["nanopore_tools"]
coverages = config["coverages"]
subsampling = config["subsampling"]
subsampled_reads_dir = config["subsampled_reads_dir"]
references_file = config["references"]
fast5s_file = config["fast5s"]
snippy_container = config["containers"]["snippy"]
samtools_container = config["containers"]["samtools"]
medaka_container = config["containers"]["medaka"]
nanopolish_container = config["containers"]["nanopolish"]


subsampled_reads = pd.read_csv(subsampled_reads_dir)
subsampled_reads = update_to_absolute_path(subsampled_reads, ["sample_path"])
references = pd.read_csv(references_file)
references = update_to_absolute_path(references, ["compressed_file", "uncompressed_file"])
fast5s_df = pd.read_csv(fast5s_file)
fast5s_df = update_to_absolute_path(fast5s_df, ["fast5_dir", "sequencing_summary"])


# ======================================================
# Main rule
# ======================================================
def get_final_files():
    final_files = []
    if illumina_tools is not None:
        final_files.extend(
            expand(f"{output_folder}/{{tool}}/illumina/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.vcf",
            tool=illumina_tools, coverage=coverages, sample=subsampled_reads["sample_id"], reference=references["reference_id"])
        )
        final_files.extend(
            expand(f"{output_folder}/{{tool}}/illumina/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.ref.fa",
            tool=illumina_tools, coverage=coverages, sample=subsampled_reads["sample_id"], reference=references["reference_id"])
        )
    if nanopore_tools is not None:
        final_files.extend(
            expand(f"{output_folder}/{{tool}}/nanopore/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.vcf",
            tool=nanopore_tools, coverage=coverages, sample=subsampled_reads["sample_id"], reference=references["reference_id"])
        )
        final_files.extend(
            expand(f"{output_folder}/{{tool}}/nanopore/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.ref.fa",
            tool=nanopore_tools, coverage=coverages, sample=subsampled_reads["sample_id"], reference=references["reference_id"])
        )
    return final_files

rule all:
    input: get_final_files()


include: "rules/snippy.smk"
include: "rules/samtools.smk"
include: "rules/medaka.smk"
include: "rules/nanopolish.smk"
