print("execute retrotransposons.py")

rule index_bam:
    input: 
        cbam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam"
    output:
        bai="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam.bai"
    benchmark: 
        "results/benchmarks/{group}__{sample_prep}.mitof_index_bam.json"
    conda: 
        "STAR.yaml"
    params: 
        n=config["max_cores"],
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    shell:
        "samtools index -b -@ {params.n} {input.cbam}"
        

rule count_reads:
    input: 
        bam= "results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
        bai="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam.bai"
    output: 
        reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
    log: 
        "results/logs/{group}__{sample_prep}_reads.txt"
    conda: 
        "retro_index.yaml"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    shell:
        "samtools view -c {input.bam} > {output.reads}"    

rule retro_count:
        """
        single threaded
        stranded matters, paired end is automatic

        """
        input: 
            bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
            bami="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam.bai",
            reference_gtf=features["reference"]["reference_gtf"],
            retro_gtf=features["reference"]["retro_gtf"],
            reads="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
        output: 
            "results/retrotransposons/{group}__{sample_prep}.cntTable"
        benchmark: 
            "results/benchmarks/retrotransposons/{group}__{sample_prep}.retro.txt"
        log: 
            "results/logs/retrotransposons/{group}__{sample_prep}.retro.txt"
        conda: 
            "retro.yaml"
        resources:
            load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        params: 
            n=config["max_cores"],
            stranded=config["transposon_stranded"]
        shell: "echo {wildcards.group}__{wildcards.sample_prep} ;\
                TEtranscripts --version ;\
                TEcount -b {input.bam}  \
               --GTF {input.reference_gtf}  \
               --TE {input.retro_gtf}  \
               --sortByPos \
               --stranded {params.stranded} \
               --project results/retrotransposons/{wildcards.group}__{wildcards.sample_prep}"
               # ;\mv {wildcards.group}__{wildcards.sample_prep}.cntTable results/retrotransposons/{wildcards.group}__{wildcards.sample_prep}.cntTable
               
## this one needs to be adjusted for rat
print("3")
if config["transposon_species"] == "rat":
    rule retrotransposons_cleanup:
            """
            this is supposed to cleanup the file
            """
            input: 
                "results/retrotransposons/{group}__{sample_prep}.cntTable"
            output: 
                "results/retrotransposons/clean/{group}__{sample_prep}.cntTable"
            shell: "awk 'NR >= 32884 && NR <= 34072' {input} > {output}"
            
    rule retrotransposons_genex:
            """
            this is supposed to cleanup the file
            """
            input: 
                "results/retrotransposons/{group}__{sample_prep}.cntTable"
            output: 
                "results/retrotransposons/genex/clean/{group}__{sample_prep}.cntTable"
            shell: "awk 'NR >= 2 && NR <= 32884' {input} > {output}"            


              

            
rule retro_tetranscripts_cleanup:
    """
    should use wc -l ce11_rmsk_TE.nu.gtf (retro_gtf)
    """
    input: 
       in1="results/retrotransposons/{group}__{sample_prep}.cntTable"
    output: 
      out1="results/retrotransposons/clean/{group}__{sample_prep}.cntTable",
      out2="results/retrotransposons/genex/clean/{group}__{sample_prep}.cntTable"
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
            
            

####################################################################################      
####################################################################################          
rule retro_summary_stats:
        """
        this 
        """
        input: 
            in1="results/retrotransposons/clean/{group}__{sample_prep}.cntTable",
            in2="results/retrotransposons/read_counts/{group}__{sample_prep}_reads.txt",
            in3="results/retrotransposons/genex/clean/{group}__{sample_prep}.cntTable",
        output: 
            "results/retrotransposons/clean/summaries/{group}__{sample_prep}.csv"
        conda: 
            "retro2.yaml"
        shell:
            "Rscript src/scripts/retro_summary_stats.R {input.in1} {output} {input.in2} {input.in3}"  
               
rule retro_cbind_filter :
        """
        this 
        """
        input: 
            in1=expand("results/retrotransposons/clean/summaries/{group}__{sample_prep}.csv",zip, group=GROUPS, sample_prep = SAMPLES_prep),
            in2=expand("results/retrotransposons/genex/clean/{group}__{sample_prep}.cntTable",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            "results/retrotransposons/final/final_clean_pooled_result.csv",
            "results/retrotransposons/final/genex_final_clean_pooled_result.csv",
        params:
            read_cutoff=config["transposon_read_cutoff"]
        conda: 
            "retro2.yaml"
        script: "../scripts/retro_cbind_filter.R"              
            
            
              
