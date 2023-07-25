library(data.table) # for fread()
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

gtf_to_filter=args[1] # "/home/kamil/bio/snakemake/results/nearest/human_refgtf_transcript_only.bed"
out_file <- args[2] # "/home/kamil/bio/snakemake/results/nearest/TEST human_refgtf_transcript_only.bed"

gtf <- fread(gtf_to_filter, sep="\t", header=F)
#head(gtf)
dim(gtf)

gtf$len <- abs(gtf$V2-gtf$V3)

# V4 is ENSEMBL ID
gtf_new=gtf %>%
  group_by(V4)%>%
  filter(len == max(len)) %>%
  ungroup()
dim(gtf_new)

print(out_file)
write.table(gtf_new[,1:10], out_file, quote=F, row.names = F, sep="\t", col.names = F)# sep="\t"
