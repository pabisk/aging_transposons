rule bwa_index_mod:
    """
    alpha
    """
    input:
        genome_fasta=features["reference"]["genome_fasta"],
    output:
        genome_fasta=features["reference"]["genome_fasta"]+".ann",
        out="results/retrotransposons_l1em/index_generated.txt"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    conda:
        "l1em.yaml"
    shell:
        "echo bwa index {input} ; \
        bwa index {input} ; \
        echo x > {output.out} "


rule run_l1em_mod:
    """
    alpha, bwa index multithreading?
    """
    input:
        bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
        bai="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam.bai",
        reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
        genome_fasta2=features["reference"]["genome_fasta"]+".ann",
        genome_fasta=features["reference"]["genome_fasta"],
    output:
        "results/retrotransposons_l1em_mod/n_{group}__{sample_prep}/L1EM-modified/full_counts.txt"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
        pwd=os.getcwd(),
        n=config["l1em_max_cores"],
        species=config["l1em_species"],
        stranded=config["l1em_stranded"],
        #l1em_path=config["l1em_path"],
    conda:
        "l1em.yaml"
    log:
        "results/logs/retrotransposons_l1em_mod/n_{group}__{sample_prep}/l1em.txt"
    shell:
        "cp -r data/L1EM-modified results/retrotransposons_l1em_mod/n_{wildcards.group}__{wildcards.sample_prep} ; \
        cd results/retrotransposons_l1em_mod/n_{wildcards.group}__{wildcards.sample_prep}/L1EM-modified ; \
        pwd ; \
        rm -d -r -f idL1reads ; \
        rm -d -r -f split_fqs ; \
        rm -d -r -f G_of_R ; \
        echo now generate fasta index: ; \
        bash generate{params.species}_L1EM_fasta_and_index.sh {params.pwd}/{input.genome_fasta} ; \
        echo now run l1em: ; \
        bash -e run_L1EM{params.species}.sh {params.pwd}/{input.bam} {params.pwd}/results/retrotransposons_l1em_mod/n_{wildcards.group}__{wildcards.sample_prep}/L1EM-modified {params.pwd}/{input.genome_fasta} {params.n} ; \
        rm -d -r -f idL1reads ; \
        rm -d -r -f split_fqs ; \
        rm -d -r -f G_of_R"

rule l1em_rbind_filter_mod:
        """
        this 
        """
        input: 
            in1=expand("results/retrotransposons_l1em_mod/n_{group}__{sample_prep}/L1EM-modified/full_counts.txt",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            "results/retrotransposons_l1em_mod/l1em_final_result.csv",
            "results/retrotransposons_l1em_mod/l1em_clean_final_result.csv",
        params:
            read_cutoff_fold=config["l1em_min_reads"],
            read_cutoff_mind_reads=config["l1em_min_fold_enrich"],
        conda: 
            "retro2.yaml"
        script: "../scripts/l1em_rbind_filter.R" 
        
        
#### experimental
