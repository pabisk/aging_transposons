library(Rsubread)
args = commandArgs(trailingOnly=TRUE)
#stranded=config["artdeco_stranded_htseq"], # 0,1,2 in featurecounts, yes, no, reverse in config
#paired=config["artdeco_type"], # T vs F in featurecounts,  SE or PE in config
#{input.gtf} {params.n} {params.stranded} {params.paired} {input.bam} 
gtf_file <- args[1] 
threads <-args[2]
strand_type <-args[3]
paired <-args[4]
out_file <-args[5]
## careful with adding new arguments, check below samples <- args[7:length(args)]
f_overlap <-as.numeric(args[6])
rt_lens <- read.csv(args[7], sep="\t", header=F)

if(paired=="PE"){paired <- TRUE}
if(paired=="SE"){paired <- FALSE}
if(strand_type=="yes"){strand_type <-1}
if(strand_type=="no"){strand_type <-0}
if(strand_type=="reverse"){strand_type <-2}

print("paired")
print(paired)

print("strand_type")
print(strand_type)

# give samples via snakemake expand
#samples <- args[7:8]
samples <- args[8:length(args)]
print("samples")
print(samples)
rez <- featureCounts(samples,
              annot.ext=gtf_file, isGTFAnnotationFile=T,
              GTF.featureType="exon", isPairedEnd=paired, nthreads=threads, strandSpecific=strand_type,
              fracOverlap=f_overlap)

## we use the original bed file from artdeco to annotate the output with RT region lengths
#rt_lens <- read.csv("/home/kamil/bio/snakemake/results/ARTDeco-master/preprocess_files/readthrough.bed", sep="\t", header=F)
rt_lens$len <- abs(rt_lens$V2-rt_lens$V3)

nurez<-data.frame(cbind(rownames(rez$counts),rownames(rez$counts),rez$counts))

rownames(nurez)<-NULL
print(head(nurez))
colnames(nurez)[1] <- "ID"
colnames(nurez)[2] <- "ID2"

rt_lens$ID2 <- rt_lens$V4
merged_df <- merge(nurez, rt_lens, by="ID2", all.x=T, all.y=F, sort=F)

nurez$ID2 <- merged_df$len

write.csv(nurez, out_file, quote=F, row.names = F)
