library("data.table")

args = commandArgs(trailingOnly=TRUE)
gtf <- args[1]
introns <- args[2]
out <- args[3]
species <- args[4]

gtf <- paste(getwd(),"/results/nearest/",species,"_refgtf.bed", sep="") #"/home/kamil/Desktop/n/snakemake_transposons2022/
print(gtf)
print("running get_intron_strandedness.R")

all_genes_bed <- fread(gtf, sep = "\t", header=F) 
all_introns_bed <- fread(introns, sep = "\t", header=F) 
head(all_genes_bed)


stranded_coord <- data.frame(name=all_genes_bed$V4, strand=all_genes_bed$V6)[!duplicated(all_genes_bed$V4),]

stranded_coord_subset <- stranded_coord[stranded_coord$name %in% all_introns_bed$V4,]
stranded_coord_subset <- stranded_coord_subset[stranded_coord_subset$name %in% all_introns_bed$V4,]
colnames(all_introns_bed) <- c("chr","start","end","name","len")

merged <- merge(x=all_introns_bed, y=stranded_coord_subset, by="name",all.x=TRUE, sort=F)
dim(merged)
dim(all_introns_bed)
dim(stranded_coord_subset)

merged <- cbind(merged$chr, merged$start, merged$end, merged$name, merged$len, merged$strand)
print(head(merged))
write.table(merged, out,
            sep = "\t", col.names = F, row.names = F, quote = FALSE)
