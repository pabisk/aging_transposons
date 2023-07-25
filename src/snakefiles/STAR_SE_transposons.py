rule STAR_index:
    input:
        genome_fasta=features["reference"]["genome_fasta"],
        reference_gtf = features["reference"]["reference_gtf"]
    output:
        directory(STAR_index_dir)
    benchmark: 
        "results/benchmarks/STAR_index.txt"
    log: 
        "results/logs/STAR_index.txt"
    conda: 
        "STAR.yaml"
    params: 
        n=config["max_cores"]
    threads: config["max_cores"]
    resources:
        load = ( lambda wildcards, attempt: min(config["star_load"] * (2 ** (attempt-1)), 100) )
    shell: 
        "mkdir -p {output}  ; \
        STAR  \
        --runThreadN {params.n}  \
        --runMode genomeGenerate  \
        --genomeDir {output}  \
        --sjdbGTFfile {input.reference_gtf}  \
        --genomeFastaFiles {input.genome_fasta} 2> {log}"


rule STAR_align_SE:
    input: 
        read1="results/trimmed_fastq/{group}__{sample_prep}.fq.gz",
        index = STAR_index_dir,
        reference_gtf = features["reference"]["reference_gtf"]
    output: 
        bam= "results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
    benchmark: 
        "results/benchmarks/{group}__{sample_prep}.align.json"
    log: 
        "results/logs/STAR/{group}__{sample_prep}.align.txt"
    conda: 
        "STAR.yaml"
    params: 
        n=config["max_cores"],
        f_multimap=config["transposon_star_outFilterMultimapNmax"],
        a_multimap=config["transposon_star_winAnchorMultimapNmax"],
        pwd=os.getcwd(),
	ul=config["star_ulimit"]
    threads: config["max_cores"]
    resources:
        load = ( lambda wildcards, attempt: min(config["star_load"] * (2 ** (attempt-1)), 100) )
    shell: "ulimit -n {params.ul} && \
		STAR \
		--runThreadN {params.n} \
		--genomeDir {input.index} \
		--readFilesIn {input.read1} \
		--outFileNamePrefix results/final_mapping/{wildcards.group}__{wildcards.sample_prep}. \
 		--readFilesCommand pigz -dc \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax {params.f_multimap}  \
		--winAnchorMultimapNmax {params.a_multimap} \
		--alignInsertionFlush Right \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 8 \
		--chimOutJunctionFormat 1 \
		--chimOutType Junctions WithinBAM \
		--chimMultimapScoreRange 3 \
		--chimScoreJunctionNonGTAG -4 \
		--chimMultimapNmax 20 \
		--chimNonchimScoreDropMin 10 \
		--alignMatesGapMax 100000 \
		--sjdbGTFfile {input.reference_gtf} \2 > {log}"
