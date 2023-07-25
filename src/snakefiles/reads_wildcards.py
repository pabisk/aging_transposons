from glob import glob
from os.path import join
from os import listdir
from os.path import isfile, join

print("executing reads_wildcards.py")
STAR_index_dir = "results/STAR_index"
reads_fq = "data/reads_fq"

FILE_ENDING = "fq.gz"
mypath=reads_fq
fq_ending=0
fastq_ending=0
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print("onlyfiles")
print(onlyfiles)

for x in onlyfiles:
    print(x)
    fq_ending=fq_ending+x.count("fq.gz")
    fastq_ending=fastq_ending+x.count("fastq.gz")

if fq_ending > fastq_ending:
    GROUPS,SAMPLES_prep, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_1.fq.gz'))
    #GROUPS2,SAMPLES_prep2, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_2.fq.gz'))
    FILE_ENDING = "fq.gz"
if fq_ending < fastq_ending:
    GROUPS,SAMPLES_prep, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_1.fastq.gz'))
    #GROUPS2,SAMPLES_prep2, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_2.fastq.gz'))
    FILE_ENDING = "fastq.gz"

if fq_ending == fastq_ending:
    print("ERROR - most likely mismatched file endings")

## below is obsolete??    
#genome_wc, = glob_wildcards(join(features["reference"]["genome_fasta"], '{genome}.genome.fa'))
#PATTERN_R1_prep = '{sample_prep}_1.fq.gz'
#PATTERN_R2_prep = '{sample_prep}_2.fq.gz'from glob import glob
from os.path import join
from os import listdir
from os.path import isfile, join

print("executing reads_wildcards.py")
STAR_index_dir = "results/STAR_index"
reads_fq = "data/reads_fq"

FILE_ENDING = "fq.gz"
mypath=reads_fq
fq_ending=0
fastq_ending=0
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print("onlyfiles")
print(onlyfiles)

for x in onlyfiles:
    print(x)
    fq_ending=fq_ending+x.count("fq.gz")
    fastq_ending=fastq_ending+x.count("fastq.gz")

if fq_ending > fastq_ending:
    GROUPS,SAMPLES_prep, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_1.fq.gz'))
    #GROUPS2,SAMPLES_prep2, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_2.fq.gz'))
    FILE_ENDING = "fq.gz"
if fq_ending < fastq_ending:
    GROUPS,SAMPLES_prep, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_1.fastq.gz'))
    #GROUPS2,SAMPLES_prep2, = glob_wildcards(join(reads_fq, '{group}__{sample_prep}_2.fastq.gz'))
    FILE_ENDING = "fastq.gz"

if fq_ending == fastq_ending:
    print("ERROR - most likely mismatched file endings")

## below is obsolete??    
#genome_wc, = glob_wildcards(join(features["reference"]["genome_fasta"], '{genome}.genome.fa'))
#PATTERN_R1_prep = '{sample_prep}_1.fq.gz'
#PATTERN_R2_prep = '{sample_prep}_2.fq.gz'
