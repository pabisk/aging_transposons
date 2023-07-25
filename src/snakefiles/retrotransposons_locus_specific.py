print("executing retrotransosons_locus_specific.py")
rule namesort:
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
    output:
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByName.out.bam",
    conda:
        "samtools.yaml"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params: 
        n=config["max_cores"]
    shell:
        "samtools sort {input.bam} -@{params.n} -n -o {output.bam}"


rule install_mod_te_trans:
    """
    verbose mode in this fork is broken, need to namesort prior to running
    CAVE: this will need to be modified when running multiple samples
    This uses namesort whereas normal tet uses coord sort
    """
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByName.out.bam",
	    reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
    output:
	    "results/retrotransposons_locus/{group}__{sample_prep}.ele.cntTable",
	    "results/retrotransposons_locus/{group}__{sample_prep}.instance.cntTable",
	    "results/retrotransposons_locus/{group}__{sample_prep}.discount.instance.cntTable",
	    "results/retrotransposons_locus/{group}__{sample_prep}.discount.ele.cntTable",
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    stranded=config["transposon_stranded"],
	    retro_locus_gtf=features["reference"]["retro_locus_gtf"],
	    retro_locus_te=features["reference"]["retro_locus_te"],
	    retro_locus_intron=features["reference"]["retro_locus_intron"],
	    retro_locus_exon_te=features["reference"]["retro_locus_exon_te"],
	    retro_locus_intergenic_te=features["reference"]["retro_locus_intergenic_te"],
    log:
        "results/logs/retrotransposons_locus/log_{group}__{sample_prep}.cntTable",
    conda:
        "retrotransposons_locus_specific.yaml"
    shell:
        "cd data/HanLabUNLV_tetoolkit/tetoolkit-master ; \
        echo $PATH ; \
        echo echo_end ; \
        echo PATH=$PATH:{params.pwd}/data/HanLabUNLV_tetoolkit/installed/bin ; \
        echo x ; \
        echo y ; \
        export PATH=$PATH:{params.pwd}/data/HanLabUNLV_tetoolkit/installed/bin ; \
        export PYTHONPATH={params.pwd}/data/HanLabUNLV_tetoolkit/installed/lib/python2.7/site-packages ; \
        echo this is the path variable: ; \
        echo $PATH ; \
        python setup.py install --prefix {params.pwd}/data/HanLabUNLV_tetoolkit/installed ; \
        TEtranscripts --version ;\
        echo x > {params.pwd}/results/toolkit_installed.txt ;\
	echo TOOLKIT INSTALLED ; \
        TEtranscripts -t {params.pwd}/{input.bam}  \
        --GTF {params.pwd}/{params.retro_locus_gtf} \
        --TE {params.pwd}/{params.retro_locus_te} \
        --intron {params.pwd}/{params.retro_locus_intron} \
        --exonTE {params.pwd}/{params.retro_locus_exon_te} \
        --intergenicTE {params.pwd}/{params.retro_locus_intergenic_te} \
        --format BAM \
        --stranded {params.stranded} \
        --project {params.pwd}/results/retrotransposons_locus/{wildcards.group}__{wildcards.sample_prep} | tee {params.pwd}/{log}"
        

if config["transposon_species"] == "human":
    rule retro_locus__cleanup2:
            """
            this is supposed to cleanup the file
            """
            input: 
                in1="results/retrotransposons_locus/{group}__{sample_prep}.instance.cntTable",
                in2="results/retrotransposons_locus/{group}__{sample_prep}.ele.cntTable",
                in3="results/retrotransposons_locus/{group}__{sample_prep}.discount.instance.cntTable",
                in4="results/retrotransposons_locus/{group}__{sample_prep}.discount.ele.cntTable",
            output: 
                out1="results/retrotransposons_locus/instance/{group}__{sample_prep}.instance.cntTable",
                out2="results/retrotransposons_locus/genex/{group}__{sample_prep}.genex.cntTable",
                out3="results/retrotransposons_locus/instance/{group}__{sample_prep}.discount.instance.cntTable",
                out4="results/retrotransposons_locus/ele/{group}__{sample_prep}.ele.cntTable",
                out5="results/retrotransposons_locus/ele/{group}__{sample_prep}.discount.ele.cntTable",
            shell:
                "awk 'NR >= 28056 && NR <= 6581396' {input.in1} > {output.out1} ; \
                awk 'NR >= 2 && NR <= 28055' {input.in2} > {output.out2} ; \
                awk 'NR >= 28056 && NR <= 6581396' {input.in3} > {output.out3} ; \
                awk 'NR >= 28056 && NR <= 30422' {input.in2} > {output.out4} ; \
                awk 'NR >= 2 && NR <= 1187' {input.in4} > {output.out5} "
       
if config["transposon_species"] == "rat":
    rule retro_locus__cleanup2:
            """
            COORDINATES NEED TO BE UPDATED!!!!
            """
            input: 
                in1="results/retrotransposons_locus/{group}__{sample_prep}.instance.cntTable",
                in2="results/retrotransposons_locus/{group}__{sample_prep}.ele.cntTable",
                in3="results/retrotransposons_locus/{group}__{sample_prep}.discount.instance.cntTable",
                in4="results/retrotransposons_locus/{group}__{sample_prep}.discount.ele.cntTable",
            output: 
                out1="results/retrotransposons_locus/instance/{group}__{sample_prep}.instance.cntTable",
                out2="results/retrotransposons_locus/genex/{group}__{sample_prep}.genex.cntTable",
                out3="results/retrotransposons_locus/instance/{group}__{sample_prep}.discount.instance.cntTable",
                out4="results/retrotransposons_locus/ele/{group}__{sample_prep}.ele.cntTable",
                out5="results/retrotransposons_locus/ele/{group}__{sample_prep}.discount.ele.cntTable",
            shell:
                "awk 'NR >= 28056 && NR <= 6581396' {input.in1} > {output.out1} ; \
                awk 'NR >= 2 && NR <= 28055' {input.in2} > {output.out2} ; \
                awk 'NR >= 28056 && NR <= 6581396' {input.in3} > {output.out3} ; \
                awk 'NR >= 28056 && NR <= 30422' {input.in2} > {output.out4} ; \
                awk 'NR >= 2 && NR <= 1187' {input.in4} > {output.out5} "

rule retro_locus_summary_stats2:
	"""
	this produces an element level summary for a single sample
	"""
	input: 
		in1="results/retrotransposons_locus/ele/{group}__{sample_prep}.ele.cntTable",
		in2="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
		in3="results/retrotransposons_locus/ele/{group}__{sample_prep}.ele.cntTable",
		in4="results/retrotransposons_locus/ele/{group}__{sample_prep}.discount.ele.cntTable"
	output: 
		"results/retrotransposons_locus/clean/summaries/{group}__{sample_prep}.csv"
	conda: 
		"retro2.yaml"
	shell:
		"Rscript src/scripts/retro_summary_stats_locus.R {input.in1} {output} {input.in2} {input.in3} {input.in4}"  

        
rule retro_locus_cbind_filter :
        """
        this 
        """
        input: 
            in1=expand("results/retrotransposons_locus/clean/summaries/{group}__{sample_prep}.csv",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in2=expand("results/retrotransposons_locus/genex/{group}__{sample_prep}.genex.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in3=expand("results/retrotransposons_locus/instance/{group}__{sample_prep}.instance.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            "results/retrotransposons_locus/final/ele_final_clean_pooled_result.csv",
            "results/retrotransposons_locus/final/genex_final_clean_pooled_result.csv",
            "results/retrotransposons_locus/final/instance_final_clean_pooled_result.csv",
        params:
            read_cutoff=config["transposon_read_cutoff"]
        conda: 
            "retro2.yaml"
        script: "../scripts/retro_cbind_filter_locus.R"  
