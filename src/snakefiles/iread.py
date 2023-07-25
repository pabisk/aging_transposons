print("executing iread.py")

rule clean_intron_bed:
    """
    """
    input: 
	    reference_gtf = features["reference"]["reference_gtf"],
	    introns="results/nearest/introns.bed",
    output:
	    introns="results/nearest/introns2.bed",
    resources:
        load = 100
    shell:
    	"""
    	chrc=$(grep -c chr {input.reference_gtf})
    	echo $chrc
    	if (( chrc > 100 )); then
    	echo gtf contains chr string, therefore we will add chr prefix to the intron.bed
    	sed 's/^/chr/' {input.introns} > {output.introns}
    	fi
    	"""
    	
rule install_and_run_ireads:
    """
        	-t is just  for FPKM calculation
        	may need to install sudo apt-get install -y libparallel-forkmanager-perl
        	multi threading comment: instead of running 1 instance with many threads we run multiple instances with 2 threads per instance
    """
    input: 
	    bam="results/final_mapping/{group}__{sample_prep}.Aligned.sortedByCoord.out.bam",
	    reference_gtf = features["reference"]["reference_gtf"],
	    introns="results/nearest/introns2.bed",
    output:
	    out1="results/retained_introns/{group}__{sample_prep}/{group}__{sample_prep}.Aligned.sortedByCoord.out.ir.txt",
	    out2=directory("results/retained_introns/{group}__{sample_prep}/")
    resources:
        load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
    params:
	    pwd=os.getcwd(),
    conda:
        "iread.yaml"
    shell:
    	"""
    	cd {params.pwd}/data/iread-master
    	export PATH=$PATH:{params.pwd}/data/iread-master
    	#export PATH=$PATH:/home/kamil/miniconda3/envs/iread/bin/python
    	export PYTHONPATH={params.pwd}/data/iread-master
    	#export PYTHONPATH=/home/kamil/miniconda3/envs/iread/bin/python
    	chmod +x count_intronic_reads.pl
    	chmod +x assess_intron.py
    	chmod +x bam2intron
    	chmod +x bedscore
    	cpan Parallel::ForkManager
    	which perl
    	export PERL5LIB=/usr/bin
    	echo $PERL5LIB
    	echo install finished, running iread.py:
    	echo python2 iread.py {params.pwd}/{input.bam} {params.pwd}/{input.introns} -o {params.pwd}/{output.out2} -t 999 --threads 2 --n_cores 2
    	python2 iread.py {params.pwd}/{input.bam} {params.pwd}/{input.introns} -o {params.pwd}/{output.out2} -t 999 --threads 2 --n_cores 2
    	"""
    	
    	
    	
rule iread_filter :
        """
        this 
        """
        input: 
            in1=expand("results/retained_introns/{group}__{sample_prep}/{group}__{sample_prep}.Aligned.sortedByCoord.out.ir.txt",zip, group=GROUPS, sample_prep = SAMPLES_prep),
        output: 
            "results/retained_introns/final_result.csv",
        params:
            read_cutoff_fold=config["l1em_min_reads"],
            read_cutoff_mind_reads=config["l1em_min_fold_enrich"],
        conda: 
            "retro2.yaml"
        script: "../scripts/iread_cbind_filter.R" 
