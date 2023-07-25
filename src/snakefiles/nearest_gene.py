print("executing telocal.py")
if config["transposon_species"] == "human" or config["transposon_species"] == "mouse":
  print("transposon_species is human")
  cleanup_string = "grep chr |"
if config["transposon_species"] == "elegans":
  print("transposon_species is elegans")
  cleanup_string = ""
  
  #### TO DO
### implement gtf2bed < /home/kamil/Desktop/n/snakemake_transposons2022/data/mouse_reference/GRCm39_rmsk_TE.gtf > test.bed
## then check: dim(tsubset)[1]/dim(res_sig_transposons)[1] on mouse dataset
## check difference between data/mouse_reference/GRCm39refGene.gtf and mm10_rmsk_TE.gtf

rule find_nearest_gene2gene:
    """
    finds the nearest gene to each gene excluding itself via ignore overlaps
    not working very well, implemented in R instead
    """
    input: 
      genes=features["reference"]["reference_gtf"],
    output:
      "results/nearest/f/find_nearest_gene2gene_completed.txt",
    resources:
      load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
      pwd=os.getcwd(),
      sp=config["transposon_species"],
      clstring=cleanup_string
    conda:
      "nearest_gene2gene.yaml"
    shell: #cat {input.genes} | grep chr | grep transcript_id | grep gene_id | convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_refgtf.bed ; \
      """
      echo step1: generate a bed with only transcripts and one with transcripts, exons, etc.
      cat {input.genes}  | {params.clstring} grep transcript_id | grep gene_id | grep -P 'transcript\t' | convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only2.bed
      echo step2: find the longest transcript with a DIY script
      Rscript src/scripts/select_longest_transcript.R {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only2.bed  {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed
      echo step 3: bedtools closest bed with itself
      mkdir -p {params.pwd}/results/nearest/gene2gene/
      echo bedtools closest -fd -t first -D ref -a {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed -b {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed
      bedtools closest -fu -t first -D ref -a {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed -b {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed > {params.pwd}/results/nearest/gene2gene/{params.sp}_gene2genedist.bed
      echo done
      rm {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only2.bed
      
      echo x > {output}
      """
      # slow as fuck, rm {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed
      # cgat gtf2gtf --method=filter --filter-method=longest-transcript {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only2.bed > {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only3.bed
      
rule find_nearest_gene:
    """
    find the nearest gene for each transposon, exclude nonstandard chromosomes and patches
    this provides transposons that are WITHIN genes (distance 0) as well as those that are upstream or downstream
    -D signed distance
    -t resolve ties
    exclude exons, only include transcripts: grep -P 'transcript\t'
    """
    input: 
      genes=features["reference"]["reference_gtf"],
      transposons=features["reference"]["retro_gtf"],
    output:
      "results/nearest/f/find_nearest_gene_completed.txt",
    resources:
      load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
      pwd=os.getcwd(),
      sp=config["transposon_species"],
      clstring=cleanup_string
    conda:
      "nearest.yaml"
    shell: #cat {input.genes} | grep chr | grep transcript_id | grep gene_id | convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_refgtf.bed ; \
      """
      echo step1: generate a bed with only transcripts and one with transcripts, exons, etc.
      cat {input.genes}  | {params.clstring} grep transcript_id | grep gene_id | grep -P 'transcript\t' | convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only.bed
      cat {input.genes}  | {params.clstring} grep transcript_id | grep gene_id | convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_refgtf.bed
      echo step2
       cat {input.transposons} | {params.clstring} convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_transposons.bed
      # cat {input.transposons} | {params.clstring} gtf2bed > {params.pwd}/results/nearest/{params.sp}_transposons.bed
      # cat {input.transposons} | {params.clstring} gtf2bed > {params.pwd}/results/nearest/{params.sp}_transposons.bed
      echo step3
      bedtools closest -t first -D b -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only.bed > {params.pwd}/results/nearest/{params.sp}_intersection.bed
      echo done
      echo x > {output}
      """

rule find_transp_in_intron:
    """
    find out which transposons are  located in introns
    -f 1 forces 100% overlap of transposons, ie transposon must be in intron
    this rule depends on find_nearest_gene because it uses the cleaned up bed file
	Software: bedtools, gtftools
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    "results/nearest/f/intersect_intron_transposon_completed.txt",
	    "results/nearest/introns.bed"
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
      "nearest.yaml"
    shell:
    	"python3 data/GTFtools_0.8.5/gtftools.py -d {params.pwd}/results/nearest/introns.bed {input.gtf} ; \
    	sed 's/chr//g' {params.pwd}/results/nearest/{params.sp}_transposons.bed > {params.pwd}/results/nearest/{params.sp}_transposons.chr.removed.bed ; \
    	bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.chr.removed.bed -b {params.pwd}/results/nearest/introns.bed > {params.pwd}/results/nearest/{params.sp}_intersection2.bed ; \
    	echo x > {output}"
    	# bedtools intersect -f 1 -wo -a rmsk2.bed -b intron_human_3877v2.bed > transpon_intersect_gene4.bed

rule generate_stranded_introns:
    """
    recover strand information from gtftools generated intron bed
    """
    input: 
      c="results/nearest/f/find_nearest_gene_completed.txt",
      introns="results/nearest/introns.bed",
      gtf=features["reference"]["reference_gtf"],
      transposons=features["reference"]["retro_gtf"],
    output:
      "results/nearest/introns_stranded.bed"
    params:
      pwd=os.getcwd(),
      sp=config["transposon_species"],
    conda:
      "nearest_r.yaml"
    shell:
      "Rscript src/scripts/get_intron_strandedness.R {params.pwd}/results/nearest/{params.sp}_refgtf.bed {input.introns} {output} {params.sp}"

rule find_transp_in_intron_stranded:
    """
    STRANDED VERSION!
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
	    introns="results/nearest/introns_stranded.bed"
    output:
	    "results/nearest/f/intersect_intron_transposon_stranded_completed.txt",
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"
    shell:
    	"sed 's/chr//g' {params.pwd}/results/nearest/{params.sp}_transposons.bed > {params.pwd}/results/nearest/{params.sp}_transposons.chr.removed.bed ; \
    	bedtools intersect -s -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.chr.removed.bed -b {params.pwd}/results/nearest/introns_stranded.bed > {params.pwd}/results/nearest/{params.sp}_stranded_intersection2.bed ; \
    	echo x > {output}"
    	# bedtools intersect -f 1 -wo -a rmsk2.bed -b intron_human_3877v2.bed > transpon_intersect_gene4.bed    	

############################################################################################################################
################################################# UTRS and Exons analysis
############################################################################################################################
rule generate_UTR_exon_beds:
    """
    
    """
    input: 
      c="results/nearest/f/find_nearest_gene_completed.txt",
      gtf=features["reference"]["reference_gtf"],
    output:
      "results/nearest/f/UTR_exon_beds_completed.txt",
    params:
      pwd=os.getcwd(),
      sp=config["transposon_species"]
    conda:
      "nearest_r.yaml"
    shell:
      "mkdir -p {params.pwd}/results/nearest/UTR/ ; \
      grep 5_UTR {params.pwd}/results/nearest/{params.sp}_refgtf.bed >  {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.temp.bed ; \
      grep $'UTR\t' {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.temp.bed >  {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.bed ; \
      rm {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.temp.bed ; \
      grep 3_UTR {params.pwd}/results/nearest/{params.sp}_refgtf.bed >  {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.temp.bed ; \
      grep $'UTR\t' {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.temp.bed >  {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.bed ; \
      rm {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.temp.bed ; \
      grep exon {params.pwd}/results/nearest/{params.sp}_refgtf.bed >  {params.pwd}/results/nearest/UTR/{params.sp}_exon_refgtf.bed ; \
      echo x > {output} "
      
## 
      
rule find_random_transp_in_UTR_exon_beds_unstranded:
    """
    find out which transposons are  located in UTR_exon_beds, UNstranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	"results/nearest/f/UTR_exon_beds_completed.txt",
    	"results/nearest/f/randomized_transposons_completed.txt",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    "results/nearest/f/UTR_exon_beds_random_intersect_completed.txt",
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"  
    shell: # bedtools intersect -f 1 -wao -a
    	"bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons_shuffled.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_rand_transp_in_5utr_intersection2.bed ; \
    	bedtools intersect -f 1 -wao  -a {params.pwd}/results/nearest/{params.sp}_transposons_shuffled.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_rand_transp_in_3utr_intersection2.bed ; \
    	bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons_shuffled.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_exon_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_rand_transp_in_exon_intersection2.bed ; \
    	echo x > {output}"
      
rule find_transp_in_UTR_exon_beds_unstranded:
    """
    find out which transposons are  located in UTR_exon_beds, UNstranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	"results/nearest/f/UTR_exon_beds_completed.txt",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    "results/nearest/f/UTR_exon_beds_intersect_completed.txt",
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"  
    shell: # bedtools intersect -f 1 -wao -a
    	"bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_5utr_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_transp_in_5utr_intersection2.bed ; \
    	bedtools intersect -f 1 -wao  -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_3utr_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_transp_in_3utr_intersection2.bed ; \
    	bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/results/nearest/UTR/{params.sp}_exon_refgtf.bed > {params.pwd}/results/nearest/UTR/{params.sp}_transp_in_exon_intersection2.bed ; \
    	echo x > {output}"
############################################################################################################################
################################################# randomize genes
############################################################################################################################
rule randomize_transposons2:
	"""
	1. create a file with random gtf positions
	2. intersect this file with the refernce gtf
	"""
	input:
		"results/nearest/f/find_nearest_gene_completed.txt",
		transposons=features["reference"]["retro_gtf"],
		gtf=features["reference"]["reference_gtf"],
	output:
		"results/nearest/f/randomized_transposons_completed.txt",
	params:
		pwd=os.getcwd(),
		sp=config["transposon_species"],
		genes=features["reference"]["genome_fasta"],
		clstring=cleanup_string
	resources:
		load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
	conda:
		"nearest.yaml"
	shell:
		"samtools faidx {params.genes} > {params.genes}.fai ; \
		awk -v OFS='\t' {{'print $1,$2'}} {params.genes}.fai > {params.genes}.fai2 ; \
		bedtools shuffle -i {input.transposons} -excl {input.gtf} -g {params.genes}.fai2 -seed 1990 | {params.clstring} convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_transposons_shuffled_excl.bed ; \
		bedtools shuffle -i {input.transposons} -g {params.genes}.fai2 -seed 1990 | {params.clstring} convert2bed --input=gtf > {params.pwd}/results/nearest/{params.sp}_transposons_shuffled.bed ; \
		bedtools closest -t first -D b -a {params.pwd}/results/nearest/{params.sp}_transposons_shuffled.bed -b {params.pwd}/results/nearest/{params.sp}_refgtf_transcript_only.bed > {params.pwd}/results/nearest/{params.sp}_shuffled_t_intersection.bed ; \
		echo x > {output}" ## bedtools closest -t first -D b -a 




############################################################################################################################
################################################# ARTDECO, readin, readthru analysis
############################################################################################################################


rule find_transp_in_readthru_unstranded:
    """
    find out which transposons are  located in readthru bed, this analysis is unstranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	rt="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    rt="results/nearest/f/rt/{intervals}_readthru_intersect_completed.txt",
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"
    shell:
    	"bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/{input.rt} > {params.pwd}/results/nearest/{params.sp}_transp_in_readthru_intersection2.bed ; \
    	echo x > {output}"
    	
rule find_transp_in_readin_unstranded:
    """
    find out which transposons are  located in readthru bed, this analysis is unstranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    ri=temp("results/nearest/f/readin_intersect_completed.txt"),
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"
    shell:
    	"bedtools intersect -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/{input.ri}  > {params.pwd}/results/nearest/{params.sp}_transp_readin_intersection2.bed ; \
    	echo x > {output}"

rule find_transp_in_readthru_stranded:
    """
    find out which transposons are  located in readthru bed, this analysis is stranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	rt="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/{intervals}_readthrough.bed",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    rt="results/nearest/f/rt/{intervals}_readthru_stranded_intersect_completed.txt",
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"
    shell:
    	"bedtools intersect -s -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/{input.rt} > {params.pwd}/results/nearest/{params.sp}_stranded_transp_in_readthru_intersection2.bed ; \
    	echo x > {output}"

rule find_transp_in_readin_stranded:
    """
    find out which transposons are  located in readthru bed, this analysis is stranded
    """
    input: 
    	"results/nearest/f/find_nearest_gene_completed.txt",
    	ri="results/ARTDeco-preprocess/ARTDeco-master/preprocess_files/read_in.bed",
	    gtf=features["reference"]["reference_gtf"],
	    transposons=features["reference"]["retro_gtf"],
    output:
	    ri=temp("results/nearest/f/readin_stranded_intersect_completed.txt"),
    resources:
        load = ( lambda wildcards, attempt: min(config["high_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
	    sp=config["transposon_species"]
    conda:
        "nearest.yaml"
    shell:
    	"bedtools intersect -s -f 1 -wao -a {params.pwd}/results/nearest/{params.sp}_transposons.bed -b {params.pwd}/{input.ri} > {params.pwd}/results/nearest/{params.sp}_stranded_transp_readin_intersection2.bed ; \
    	echo x > {output}"
