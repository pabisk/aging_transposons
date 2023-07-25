#############################
#shell.prefix("set -euo pipefail;")
shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

## INCLUDE START
snakefiles = "src/snakefiles/"
include: snakefiles + "snakefile_plus.py"
include: snakefiles + "folders.py"
include: snakefiles + "reads_wildcards.py"
if config["pre_process_pair_or_single_end"] == "paired":
    if config["star_fix_shortread_issues"] == "yes": 
    	print("PE + short base mode")
        include: snakefiles + "STAR_tweaked_short_reads.py"
    else:
        include: snakefiles + "STAR.py"
        print("PE + normal")
    include: snakefiles + "pre_process.py"
else:
    include: snakefiles + "pre_process_SE.py"
    include: snakefiles + "STAR_SE_transposons.py"
    print("SE + normal")
    
include: snakefiles + "retrotransposons.py"
include: snakefiles + "retrotransposons_locus_specific.py"
include: snakefiles + "l1em.py"
include: snakefiles + "telocal.py"
include: snakefiles + "nearest_gene.py"
include: snakefiles + "iread.py"
include: snakefiles + "artdeco.py"
include: snakefiles + "texp_test.py"

## INCLUDE END

rule all:
    input:
    	#### STAR only
    	expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        ####uncomment for retro locus specific analysis with tettranscripts-modified
        #"results/retrotransposons_locus/final/ele_final_clean_pooled_result.csv",
        ####uncomment for retro analysis with DEFAULT tettranscripts
        "results/retrotransposons/final/final_clean_pooled_result.csv",
        ####uncomment for l1em analysis
        #"results/retrotransposons_l1em/l1em_final_result.csv",
        #### uncomment for telocal
        "results/telocal/final/te_final_clean_pooled_result.csv",
        #### iread to check for retained introns
        #expand("results/retained_introns/{group}__{sample_prep}/{group}__{sample_prep}.Aligned.sortedByCoord.out.ir.txt",zip, group=GROUPS, sample_prep = SAMPLES_prep),
         #"results/retained_introns/final_result.csv",
         #### artdeco default and artdeco featurecounts version
         #"results/ARTDeco-master/quantification/read_in.raw.txt", # outdated
         #"results/ARTDeco-featurecounts/read_in.raw.txt", # outdated
         #"results/artdeco_finished.txt", # this is for classic artdeco
         #expand("results/ARTDeco-featurecounts/{intervals}_readthrough.raw.txt", intervals=config["artdeco_readthru_intervals"]),
         #expand("results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_{intervals}_readthrough.gtf", intervals=config["artdeco_readthru_intervals"]),
         #### simple implementation of l1em
         # expand("results/l1em-like/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
         # "results/l1em-like/final/te_final_clean_pooled_result.csv",
         #### UTR diff expression
         # expand("results/utr_diy/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
         #### find nearest genes
         "results/nearest/f/find_nearest_gene_completed.txt",
         "results/nearest/f/intersect_intron_transposon_completed.txt",
         "results/nearest/introns_stranded.bed",
         "results/nearest/f/intersect_intron_transposon_stranded_completed.txt",
         "results/nearest/f/randomized_transposons_completed.txt",
         expand("results/nearest/f/rt/{intervals}_readthru_stranded_intersect_completed.txt", intervals=config["artdeco_readthru_intervals"]),
         expand("results/nearest/f/rt/{intervals}_readthru_intersect_completed.txt", intervals=config["artdeco_readthru_intervals"]),
         "results/nearest/f/readin_intersect_completed.txt",
         "results/nearest/f/UTR_exon_beds_intersect_completed.txt",
         "results/nearest/f/UTR_exon_beds_random_intersect_completed.txt",
         "results/nearest/f/readin_stranded_intersect_completed.txt",
         "results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",

