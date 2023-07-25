rule genome_index:
    input:
        genome_fasta=features["reference"]["genome_fasta"]
    output:
        f_index=features["reference"]["genome_fasta"] + ".fai"
    conda: 
        "samtools.yaml"
    log: 
        "results/logs/fasta_index.txt"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    shell:
        "samtools faidx {input.genome_fasta}"

rule gen_chrom_sizes:
    input:
        genome_fasta=features["reference"]["genome_fasta"]+".fai",
    output:
        genome_fasta=features["reference"]["genome_fasta"]+".chrlen",
        out="results/artdeco_tempfiles/chrom_len_generated.txt"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    shell:
        "cut -f1,2 {input.genome_fasta} > {output.genome_fasta} ; \
        echo x > {output.out} "
        
rule clean_gtf:
    """
    rm non transcript ID features from gtf so it is compatible with artdeco
    """
    input:
        genome_gtf=features["reference"]["reference_gtf"]
    output:
        genome_gtf=features["reference"]["reference_gtf"]+".clean",
        out="results/artdeco_tempfiles/clean_gtf_generated.txt"
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    shell:
        "bash src/scripts/clean_gtf.sh {input.genome_gtf} {output.genome_gtf}  ; \
        echo x > {output.out}"

rule bed_to_gtf_ri:
  """
  use bed2gtf.py to reformat the artdeco outputted bed file into a format (almost) suitable for featurecounts and htseq
  """
  input:
    ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
  output:
    ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/readin.gtf",
  resources:
    load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
  params:
    pwd=os.getcwd(),
  shell:
    "echo cp {input.ri} {input.ri}.v2.bed ; \
    cp {input.ri} {input.ri}.v2.bed ; \
    python2 data/bed2gtf-master/bed2gtf.py -i {input.ri}.v2.bed -o {output.ri} ; \
    rm {input.ri}.v2.bed"
    
rule bed_to_gtf_rt:
  """
  use bed2gtf.py to reformat the artdeco outputted bed file into a format (almost) suitable for featurecounts and htseq
  """
  input:
    rt="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed",
  output:
    rt="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.gtf",
  resources:
    load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
  params:
    pwd=os.getcwd(),
  shell:
    "cp {input.rt} {input.rt}.v2.bed ; \
    echo cp {input.rt} {input.rt}.v2.bed ; \
    ls results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/ ; \
    python2 data/bed2gtf-master/bed2gtf.py -i {input.rt}.v2.bed -o {output.rt} ; \
    rm {input.rt}.v2.bed"

rule rename_gtf_ri:
  """
  rename the gtf so it is fully compatible with htseq and featurecounts, e.g. hack it to use exon attribute
  """
  input:
    gtf2="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/readin.gtf",
    bed2="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
  output:
    "results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_readin.gtf",
  resources:
    load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
  params:
    pwd=os.getcwd(),
  shell:
    "Rscript {params.pwd}/src/scripts/rename_gtf_for_artdeco.R {input.gtf2} {input.bed2}"

rule rename_gtf_rt:
  """
  rename the gtf so it is fully compatible with htseq and featurecounts, e.g. hack it to use exon attribute
  """
  input:
    gtf1="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.gtf",
    bed1="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed",
  output:
    "results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_{intervals}_readthrough.gtf",
  resources:
    load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
  params:
    pwd=os.getcwd(),
  shell:
    "Rscript {params.pwd}/src/scripts/rename_gtf_for_artdeco.R {input.gtf1} {input.bed1}"
    
rule artdeco_htseq:
  """
  this will run htseq with the bed/gtf files that artdeco generated, a hybrid mode so to say. VERY slow. artdeco takes 5h, this takes 10h for just readthru and
  featurecounts takes 20min!
  """
  input:
    gtf="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_readthrough.gtf",
    bam=expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),
  output:
    "results/ARTDeco-htseq/results/htseq_readthru.tsv",
  resources:
    load = 100
  params:
    n=config["l1em_max_cores"],
    stranded=config["artdeco_stranded_htseq"],
  conda:
    "htseq.yaml"
  shell:
    "htseq-count --stranded={params.stranded} --order=pos -n {params.n} {input.bam} {input.gtf} --counts_output={output}"
    
rule artdeco_featurecounts_ri:
  """
  this will run featurecounts with the bed/gtf files that artdeco generated, a hybrid mode so to say
  """
  input:
    gtf2="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_readin.gtf",
    bam=expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),
    ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
  output:
    ri="results/ARTDeco-featurecounts/read_in.raw.txt"
  resources:
    load = 100
  params:
    n=config["l1em_max_cores"],
    stranded=config["artdeco_stranded_htseq"], # 0,1,2 in featurecounts, yes, no, reverse in config
    paired=config["artdeco_type"], # T vs F in featurecounts,  SE or PE in config
    fracOverlap=config["featurecounts_fracOverlap"], # T vs F in featurecounts,  SE or PE in config
    len_data="results/ARTDeco-master/preprocess_files/readthrough.bed"
  log:
    l2="results/logs/artdeco_featurecounts/log_readin_artdeco_featurecounts.txt",
  conda:
    "htseq.yaml"
  shell:
    "Rscript src/scripts/featurecounts.R {input.gtf2} {params.n} {params.stranded} {params.paired} {output.ri} {params.fracOverlap} {input.ri} {input.bam}  | tee {log.l2} "
 
rule artdeco_featurecounts_rt:
  """
  this will run featurecounts with the bed/gtf files that artdeco generated, a hybrid mode so to say
  """
  input:
    gtf="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/renamed_{intervals}_readthrough.gtf",
    rt="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed"
  output:
    rt="results/ARTDeco-featurecounts/{intervals}_readthrough.raw.txt"
  resources:
    load = 100
  params:
    bam=expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),    n=config["l1em_max_cores"],
    stranded=config["artdeco_stranded_htseq"], # 0,1,2 in featurecounts, yes, no, reverse in config
    paired=config["artdeco_type"], # T vs F in featurecounts,  SE or PE in config
    fracOverlap=config["featurecounts_fracOverlap"], # T vs F in featurecounts,  SE or PE in config
    len_data="results/ARTDeco-master/preprocess_files/readthrough.bed"
  log:
    l1="results/logs/artdeco_featurecounts/{intervals}_log_readthru_artdeco_featurecounts.txt",
  conda:
    "htseq.yaml"
  shell:
    """echo {input.gtf}
    Rscript src/scripts/featurecounts.R {input.gtf} {params.n} {params.stranded} {params.paired} {output.rt} {params.fracOverlap} {input.rt} {params.bam}  | tee {log.l1}
    cp {input.rt} {output.rt}.coords.bed
    echo x > results/ARTDeco-featurecounts/readthrough.raw.txt """
    ## 

rule run_artdeco_merge:  
    input:
        expand("results/ARTDeco_{intervals}/ARTDeco-master/quantification/read_in.raw.txt", intervals=config["artdeco_readthru_intervals"]),
    output:
        "results/artdeco_finished.txt"
    shell:
        "echo x > results/artdeco_finished.txt"

rule run_artdeco:
    """
    modern libraries like trueseq are reverse stranded, this analysis wants total RNA and not poly A. NEW: added -gene-types protein_coding lincRNA
    """
    input:
        bam=expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        genome_fasta=features["reference"]["genome_fasta"],
        genome_fasta_len=features["reference"]["genome_fasta"]+".chrlen",
        genome_gtf=features["reference"]["reference_gtf"]+".clean",
    output:
        expand("results/ARTDeco_{intervals}/ARTDeco-master/quantification/read_in.raw.txt", intervals=config["artdeco_readthru_intervals_default"]),
        expand("results/ARTDeco_{intervals}/ARTDeco-master/quantification/readthrough.raw.txt", intervals=config["artdeco_readthru_intervals_default"]),
    resources:
        load = 100
    params:
        pwd=os.getcwd(),
        n=config["l1em_max_cores"],
        seqtype2=config["artdeco_type"],
        stranded=config["artdeco_stranded"],
        dist=config["artdeco_readthru_dist"],
        readin_dist=config["artdeco_readin_dist"],
        igdist=config["intergenic_max_len"],
        intervals=config["artdeco_readthru_intervals_default"],
    conda:
        "artdeco.yaml"
    log:
        expand("results/logs/artdeco/log_{intervals}_artdeco.txt", intervals=config["artdeco_readthru_intervals"]),
    shell:
        """
        echo running ARTDeco in a loop expanding the intervals
        
        IN="{params.intervals}"
        set -- "$IN" 
        IFS=" "; declare -a rt_windows=($*) 
        
        for i in "${{rt_windows[@]}}"
        do
        
        echo i
        echo $i
        echo rt windows all
        echo ${{rt_windows}}
        echo rt windows "${{rt_windows[i]}}"
        echo "${{rt_windows[i]}}"
        
        rtdist=$(echo "$i" | cut -d "-" -f 1 )
        rtdist2=$(echo "$i" | cut -d "-" -f 2 )
        
        echo rtdist and rtdist2
        echo $rtdist
        echo $rtdist2
        
        rtdist=$(($rtdist*1000))
        rtdist2=$(($rtdist2*1000))
        igdist=$(($rtdist2 - rtdist))
        
        echo ${{rt_windows[i]}}
        echo ${{rt_windows[i]}}
        
        echo cp -r {params.pwd}/data/ARTDeco-master {params.pwd}/results/ARTDeco_$i
        cp -r {params.pwd}/data/ARTDeco-master {params.pwd}/results/ARTDeco_$i
        pwd
        cd {params.pwd}/results/ARTDeco_$i/ARTDeco-master
        pwd
        python setup.py install
	
	      echo ARTDeco -mode readthrough -home-dir {params.pwd}/results/ARTDeco_$i/ARTDeco-master \
	      -bam-files-dir {params.pwd}/results/final_mapping/ \
	      -gtf-file {params.pwd}/{input.genome_gtf} -cpu {params.n} -chrom-sizes-file {params.pwd}/{input.genome_fasta_len} \
	      -layout {params.seqtype2} -stranded {params.stranded} \
	      -orientation Forward  -gene-types protein_coding lincRNA -readthrough-dist ${{rtdist}} -read-in-dist {params.readin_dist} \
	      -intergenic-max-len ${{igdist}}
	
	      ARTDeco -mode readthrough -home-dir {params.pwd}/results/ARTDeco_$i/ARTDeco-master \
	      -bam-files-dir {params.pwd}/results/final_mapping/ \
	      -gtf-file {params.pwd}/{input.genome_gtf} -cpu {params.n} -chrom-sizes-file {params.pwd}/{input.genome_fasta_len} \
	      -layout {params.seqtype2} -stranded {params.stranded} \
	      -orientation Forward  -gene-types protein_coding lincRNA -readthrough-dist ${{rtdist}} -read-in-dist {params.readin_dist} \
	      -intergenic-max-len ${{igdist}}
	      
	      cd {params.pwd}/results/ARTDeco_$i/ARTDeco-master/preprocess_files
	      GLOBIGNORE=readthrough.bed:read_in.bed
	      rm -v *.bed
	      unset GLOBIGNORE
	      
	      done
        """
        
rule run_artdeco_preprocess:
    """
    artdeco needs at least a mock bam to run and generate the intermediate files we care about, test bam will be a downscaled version of input bam
    we loop through config["artdeco_readthru_intervals"], string splitting the X-Y into -readthrough-dist and-intergenic-max-len for artdeco
    then this should generate N files equal to length of config["artdeco_readthru_intervals"]
    """
    input:
        genome_fasta=features["reference"]["genome_fasta"],
        genome_fasta_len=features["reference"]["genome_fasta"]+".chrlen",
        genome_gtf=features["reference"]["reference_gtf"]+".clean",
        dir_in=directory(STAR_index_dir),
        in1=expand("results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",zip, group=GROUPS, sample_prep = SAMPLES_prep),
    output:
        rt=expand("results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed", intervals=config["artdeco_readthru_intervals"]),
        ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
    resources:
        load = ( lambda wildcards, attempt: min(config["mid_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
        pwd=os.getcwd(),
        n=config["l1em_max_cores"],
        dist=config["artdeco_readthru_dist"],
        readin_dist=config["artdeco_readin_dist"],
        igdist=config["intergenic_max_len"],
        seqtype2=config["artdeco_type"],
        stranded=config["artdeco_stranded"],
        intervals=config["artdeco_readthru_intervals"],
    conda:
        "artdeco.yaml"
    shell:
        """
        files=( {params.pwd}/results/final_mapping/*.out.bam )
        echo "${{files[0]}}"
        samtools view -s 0.001 -b "${{files[0]}}" > {params.pwd}/data/test-bam-artdeco/test.bam
        cp -r data/ARTDeco-master results/ARTDeco-preprocess
        cd results/ARTDeco-preprocess/ARTDeco-master
        echo {params.intervals}
        python setup.py install
        echo running ARTDeco in a loop expanding the intervals
        
        IN="{params.intervals}"
        set -- "$IN" 
        IFS=" "; declare -a rt_windows=($*) 
        
        for i in "${{rt_windows[@]}}"
        do
        
        rtdist=$(echo "$i" | cut -d "-" -f 1 )
        rtdist2=$(echo "$i" | cut -d "-" -f 2 )
        
        rtdist=$(($rtdist*1000))
        rtdist2=$(($rtdist2*1000))
        igdist=$(($rtdist2 - rtdist))
	
	      echo ARTDeco -mode readthrough -home-dir {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master \
	      -bam-files-dir {params.pwd}/data/test-bam-artdeco \
	      -gtf-file {params.pwd}/{input.genome_gtf} -cpu 2 -chrom-sizes-file {params.pwd}/{input.genome_fasta_len} \
	      -layout {params.seqtype2} -stranded {params.stranded} \
	      -orientation Forward  -gene-types protein_coding lincRNA -readthrough-dist ${{rtdist}} -read-in-dist {params.readin_dist} \
	      -intergenic-max-len ${{igdist}}
	
	      ARTDeco -mode readthrough -home-dir {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master \
	      -bam-files-dir {params.pwd}/data/test-bam-artdeco \
	      -gtf-file {params.pwd}/{input.genome_gtf} -cpu 2 -chrom-sizes-file {params.pwd}/{input.genome_fasta_len} \
	      -layout {params.seqtype2} -stranded {params.stranded} \
	      -orientation Forward  -gene-types protein_coding lincRNA -readthrough-dist ${{rtdist}} -read-in-dist {params.readin_dist} \
	      -intergenic-max-len ${{igdist}}
	      
	      sleep 1
	      
	      echo mv {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/readthrough.bed {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/"$i"_readthrough.bed
	      head {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/readthrough.bed
	      mv {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/readthrough.bed {params.pwd}/results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/"$i"_readthrough.bed
	      done
        """
        
