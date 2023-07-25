print("all ok")
if config["pre_process_down_sample_rate"] < 1:
    rule downsample_se:
        input:
            read1="data/reads_fq/{group}__{sample_prep}_1."+FILE_ENDING,
        output:
            read1="results/sampled/{group}__{sample_prep}.fq.gz",
        params:
            down_sample=config["pre_process_down_sample_rate"],
            space_save=config["pre_process_disk_space_experimental"],
            end=config["pre_process_pair_or_single_end"],
            pwd=os.getcwd(),
        conda: 
            "pre_process.yaml"
        resources:
            load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        shell:
            "reformat.sh in1={input.read1} \
            out1={output.read1} \
            samplerate={params.down_sample}  ; \
            bash src/scripts/remove_input.sh {params.space_save} {params.pwd}/{input} {params.pwd}/{input} {params.end}"
## sample a fraction of fastq reads {params.down_sample}; for test runs; set to 1=100% for normal run


if config["pre_process_down_sample_rate"] == 1:
    rule fastp_se:
        input:
            read1="data/reads_fq/{group}__{sample_prep}_1."+FILE_ENDING,
        output:
            read1=temp("results/trimmed_fastq/{group}__{sample_prep}.fq.gz"),
        params: 
            n=config["max_cores"],
            phred=config["phred_score"],
            space_save=config["pre_process_disk_space_experimental"],
            pwd=os.getcwd(),
        conda: 
            "pre_process.yaml"
        log: 
            json="results/logs/fastp/{group}__{sample_prep}.fastp.json",
            htm="results/logs/fastp/{group}__{sample_prep}.fastp.html",
            general_log="results/logs/fastp/{group}__{sample_prep}.fastp.txt"
        resources:
            load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        shell:
            "fastp -i {input.read1} \
            -o {output.read1} \
            -j {log.json} \
            -h {log.htm} \
            /qualified_quality_phred {params.phred} \
            /correction \
            /verbose \
            /thread {params.n} 2> {log.general_log} ; \
            bash src/scripts/remove_input.sh {params.space_save} {params.pwd}/{input}"
## the idea is that we skip downsampling at a rate of 1 which is similar to making a copy and just takes up disk space
## performance 30M read RNAseq (3gb) files only take about 10gb RAM and 10 cores when 2 are executed in parallel

if config["pre_process_down_sample_rate"] < 1:
    rule fastp_se:
        input:
            read1="results/sampled/{group}__{sample_prep}.fq.gz",
        output:
            read1=temp("results/trimmed_fastq/{group}__{sample_prep}.fq.gz"),
        params: 
            n=config["max_cores"],
            phred=config["phred_score"],
        conda: 
            "pre_process.yaml"
        log: 
            json="results/logs/fastp/{group}__{sample_prep}.fastp.json",
            htm="results/logs/fastp/{group}__{sample_prep}.fastp.html",
            general_log="results/logs/fastp/{group}__{sample_prep}.fastp.txt"
        resources:
            load = ( lambda wildcards, attempt: min(config["low_intensity_rules_load"] * (2 ** (attempt-1)), 100) )
        shell:
            "fastp -i {input.read1} \
            -o {output.read1} \
            -j {log.json} \
            -h {log.htm} \
            /qualified_quality_phred {params.phred} \
            /correction \
            /verbose \
            /thread {params.n} 2> {log.general_log}"
## /low_complexity_filter should be DSABLED (default=off)
## /correction allow repair of paired end reads if one mate is much higher quality
## /qualified_quality_phred some phred quality cutoff ..
## adaptor trimming is on per default but I have not validated this in detail
