# aging_transposons  

In order to run the snakemake pipeline in snakemake_github/:  
data/reads_fq = put fastq files here  
data/human_reference = put annotation files in the respective folder  
Snakefile = uncomment the outputs you want  
features_homo.yaml = rename the species specific features file you will be using to just "features.yaml"  

Requires: conda, mamba, snakemake

Run with the following command, e.g.:  
time snakemake --use-conda --reason -j 8 --resources load=100 --restart-times 0 --conda-frontend mamba --keep-going  

If there is an initial error message, simply re-try with the same command.
