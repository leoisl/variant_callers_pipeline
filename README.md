# Variant callers pipeline

A snakemake pipeline for the pandora analysis pipeline that runs all the other variant callers pipelines.
Right now, it runs snippy, samtools/bcftools, medaka and nanopolish.

The version used in the pandora paper review has tag `pandora_paper_update_31_03_2021`.

The version used in the pandora paper submission has tag `pandora_paper_tag1`.

# Running

## Requirements

### Dependencies
* python 3.6+;
* singularity 2.4.1+;

### Setting up virtualenv
`./setup.sh`

## Running on the sample example:

1. Download sample data (TODO: add link);
2. `unzip sample_data.zip`
3. `source venv/bin/activate`
4. `bash scripts/run_pipeline_local.sh -j8`

## Running on the paper data:

1. `git checkout pandora_paper_tag1`
2. `source venv/bin/activate`

### If you want to run local:

3. `bash scripts/run_pipeline_local.sh -j <NB_OF_THREADS> --configfile config.pandora_paper_tag1.yaml`

### If you want to run on an LSF cluster:

3. `bash scripts/submit_lsf.sh --configfile config.pandora_paper_tag1.yaml`


### Best way to run on EBI Yoda (LSF cluster + heavy big mem node):

3.  `bash scripts/submit_lsf_nanopolish_local_bigmem_node.sh --configfile config.pandora_paper_tag1.yaml --singularity-prefix /nfs/leia/singularity/leandro/`

# Troubleshooting

If you get this error:
```
Building DAG of jobs...
Pulling singularity image docker://leandroishilima/subsampler:pandora_paper_tag1.
WorkflowError:
Failed to pull singularity image from docker://leandroishilima/subsampler:pandora_paper_tag1:
WARNING: pull for Docker Hub is not guaranteed to produce the
WARNING: same image on repeated pull. Use Singularity Registry
WARNING: (shub://) to pull exactly equivalent images.
ERROR: Image file exists, not overwriting.

  File "/hps/nobackup2/iqbal/leandro/pandora_paper_tag1/subsampler/venv/lib/python3.7/site-packages/snakemake/deployment/singularity.py", line 88, in pull
```

pass to the running script the default location where singularity images are store.
For example, in the EBI cluster, it would be `--singularity-prefix /nfs/leia/singularity/leandro/`.

