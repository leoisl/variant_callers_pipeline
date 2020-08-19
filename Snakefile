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

def get_illumina_reads(samples, sample_name, first_or_second, subsampling, coverage):
    assert first_or_second in [1, 2]
    assert sample_name in samples.sample_id.to_list()
    sample_path = samples[samples.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.{coverage}x.{subsampling}.illumina.{first_or_second}.fastq"

def get_nanopore_reads(samples, sample_name):
    assert sample_name in samples.sample_id.to_list()
    sample_path = samples[samples.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.nanopore.fastq.gz"

def get_uncompressed_reference(references, reference_id):
    assert reference_id in references.reference_id.to_list()
    uncompressed_reference_path = references[references.reference_id == reference_id]["uncompressed_file"].tolist()[0]
    return uncompressed_reference_path



# ======================================================
# Config files and vars
# ======================================================
configfile: "config.yaml"
output_folder = config['output_folder']
illumina_tools = config["illumina_tools"]
coverages = config["coverages"]
subsampling = config["subsampling"]
samples_file = config["samples"]
references_file = config["references"]
snippy_container = config["containers"]["snippy"]
samtools_container = config["containers"]["samtools"]


samples = pd.read_csv(samples_file)
samples = update_to_absolute_path(samples, ["sample_path"])
references = pd.read_csv(references_file)
references = update_to_absolute_path(references, ["compressed_file", "uncompressed_file"])


# ======================================================
# Main rule
# ======================================================
rule all:
    input:
        expand(f"{output_folder}/{{tool}}/illumina/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.vcf",
               tool=illumina_tools, coverage=coverages, sample=samples["sample_id"], reference=references["reference_id"]),
        expand(f"{output_folder}/{{tool}}/illumina/{{coverage}}x/{subsampling}/{{sample}}/{{tool}}_{{sample}}_AND_{{reference}}.ref.fa",
               tool=illumina_tools, coverage=coverages, sample=samples["sample_id"], reference=references["reference_id"]),


include: "rules/snippy.smk"
include: "rules/samtools.smk"