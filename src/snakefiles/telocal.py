print("executing telocal.py")
rule install_and_run_te_local:
    """
    BUGS lots of warning messages, try it with python 3. should be fixed
    """
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
	    reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
    output:
	    "results/telocal/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable",
    resources:
        load = ( lambda wildcards, attempt: min(config["mid_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    stranded=config["transposon_stranded"],
	    locus_gtf=features["reference"]["telocus_gtf"],
	    reference_gtf=features["reference"]["reference_gtf"],
    log:
        "results/logs/telocal/log_{group}__{sample_prep}.txt",
    conda:
        "retrotransposons_locus_specific2.yaml"
    shell:
        "cp -r -d data/TElocal-master results/telocal/{wildcards.group}__{wildcards.sample_prep}/ ; \
        cd results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master ; \
		export PATH=$PATH:{params.pwd}/results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        export PYTHONPATH={params.pwd}/results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages/ ; \
        mkdir -p {params.pwd}/results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/bin ; \
        mkdir -p {params.pwd}/results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/lib/python2.7/site-packages ; \
        echo test  ; \
        python setup.py install --prefix {params.pwd}/results/telocal/{wildcards.group}__{wildcards.sample_prep}/TElocal-master/ ; \
        echo test 3 ; \
        TElocal  --version ; \
        TElocal -b {params.pwd}/{input.bam} \
        --GTF {params.pwd}/{params.reference_gtf} \
        --stranded {params.stranded} \
        --sortByPos \
        --project {wildcards.group}__{wildcards.sample_prep} \
        --TE {params.pwd}/{params.locus_gtf} 2>&1 | tee {params.pwd}/{log}"
        
                
if config["transposon_species"] == "rat":
    rule retro_locus_cleanup_telocal:
            """
            COORDS STILL WRONG!!!!
            """
            input: 
                in1="results/telocal/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable"
            output: 
                out1="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_transposons.cntTable",
                out2="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_genex.cntTable"
            shell:
                "awk 'NR >= 17496 && NR <= 6581396' {input.in1} > {output.out1} ; \
                awk 'NR >= 2 && NR <= 3817124' {input.in1} > {output.out2}"
                

if config["transposon_species"] == "elegans" or config["transposon_species"] == "mouse" or config["transposon_species"] == "human":
    rule retro_locus_cleanup_telocal:
            """
            should use wc -l ce11_rmsk_TE.nu.gtf (retro_gtf)
            """
            input: 
                in1="results/telocal/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}.cntTable"
            output: 
                out1="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_transposons.cntTable",
                out2="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_genex.cntTable"
            params:
              reference_gtf=features["reference"]["reference_gtf"]
            shell:
                """
                wc -l {params.reference_gtf}
                wcl=$(grep -c gene$'\t' {params.reference_gtf})
                echo $wcl
                ((wcl=wcl+2))
                echo $wcl
                wcl_transposon_file_len=$(wc -l < {input.in1}) ; \
                awk -v awkvar="$wcl" -v awkvar2="$wcl_transposon_file_len" 'NR >= awkvar && NR <= awkvar2' {input.in1} > {output.out1}
                awk -v awkvar="$wcl" 'NR >= 2 && NR < awkvar' {input.in1} > {output.out2}
                """


rule retro_locus_summary_stats_telocal:
	"""
	this produces an element level summary for a single sample
	"""
	input: 
		in1="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_transposons.cntTable",
		in2="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
		#in3="results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_transposons.cntTable",
		#in4="results/retrotransposons_locus/ele/{group}__{sample_prep}.discount.ele.cntTable"
	output: 
		"results/telocal/clean/summaries/{group}__{sample_prep}.csv"
	resources:
		load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
	conda: 
		"retro2.yaml"
	shell:
		"Rscript src/scripts/retro_summary_stats_locus2.R {input.in1} {output} {input.in2}"  

rule retro_locus_cbind_filter_telocal:
        """
        this 
        """
        input: 
            in1=expand("results/telocal/clean/summaries/{group}__{sample_prep}.csv",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in2=expand("results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_genex.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in3=expand("results/telocal/clean/{group}__{sample_prep}/TElocal-master/{group}__{sample_prep}_transposons.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            "results/telocal/final/summary_te_final_clean_pooled_result.csv", 
            "results/telocal/final/genex_final_clean_pooled_result.csv",
            "results/telocal/final/te_final_clean_pooled_result.csv",
        params:
            read_cutoff=config["transposon_read_cutoff"],
            n=config["max_cores"]
        resources:
        	load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        conda:
        	"retro2.yaml"
        script: "../scripts/retro_cbind_filter_locus.R"
