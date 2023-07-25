print("executing l1em_diy.py")
rule l1em_diy:
    """
    BUGS lots of warning messages, try it with python 3. should be fixed
    """
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
	    reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
    output:
	    "results/l1em-like/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    stranded=config["transposon_stranded"],
	    locus_gtf=features["reference"]["telocus_gtf_l1em_diy"],
	    decoy_gtf=features["reference"]["decoy_gtf"],
    log:
        "results/logs/l1em-like/log_{group}__{sample_prep}.txt",
    conda:
        "retrotransposons_locus_specific2.yaml"
    shell:
        "cp -r -d data/TElocal-master results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/ ; \
        cd results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master ; \
		export PATH=$PATH:{params.pwd}/results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        export PYTHONPATH={params.pwd}/results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages/ ; \
        mkdir -p {params.pwd}/results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        mkdir -p {params.pwd}/results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages ; \
        echo test  ; \
        python setup.py install --prefix {params.pwd}/results/l1em-like/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/ ; \
        echo test 3 ; \
        TElocal  --version ; \
        TElocal -b {params.pwd}/{input.bam} \
         --GTF {params.pwd}/{params.decoy_gtf} \
        --stranded {params.stranded} \
        --sortByPos \
        --project {wildcards.group}__{wildcards.sample_prep} \
        --TE {params.pwd}/{params.locus_gtf} 2>&1 | tee {params.pwd}/{log}"
        ##


rule l1em_diy_cbind_filter_telocal:
        """
        this should include all reads for downstream calculations
        """
        input: 
            #in2=expand("results/l1em-like/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in3=expand("results/l1em-like/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            #"results/l1em-like/final/genex_final_clean_pooled_result.csv",
            "results/l1em-like/final/te_final_clean_pooled_result.csv",
        params:
            read_cutoff=0, #config["transposon_read_cutoff"]
            n=config["max_cores"]
        resources:
        	load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        conda:
        	"retro2.yaml"
        script: "../scripts/l1em_like_retro_cbind_filter_locus.R"

##############################################################
##### UTR analysis
##############################################################
rule utr_diy:
    """
    in gtf may need transcript
    """
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
	    reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
    output:
	    "results/utr_diy/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    stranded=config["transposon_stranded"],
	    locus_gtf="data/human_reference/lines_decoy.gtf.locInd",
	    #decoy_gtf=features["reference"]["decoy_gtf"],
	    decoy_gtf="data/human_reference/test3.gtf",
    log:
        "results/logs/utr_diy/log_{group}__{sample_prep}.txt",
    conda:
        "retrotransposons_locus_specific2.yaml"
    shell:
        "cp -r -d data/TElocal-master results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/ ; \
        cd results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master ; \
		export PATH=$PATH:{params.pwd}/results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        export PYTHONPATH={params.pwd}/results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages/ ; \
        mkdir -p {params.pwd}/results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        mkdir -p {params.pwd}/results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages ; \
        echo test  ; \
        python setup.py install --prefix {params.pwd}/results/utr_diy/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/ ; \
        echo test 3 ; \
        TElocal  --version ; \
        TElocal -b {params.pwd}/{input.bam} \
         --GTF {params.pwd}/{params.decoy_gtf} \
        --stranded {params.stranded} \
        --sortByPos \
        --project {wildcards.group}__{wildcards.sample_prep} \
        --TE {params.pwd}/{params.locus_gtf} 2>&1 | tee {params.pwd}/{log}"
        ##
