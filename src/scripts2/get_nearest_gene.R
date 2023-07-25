## find the nearest gene up and downstream
library("data.table")
gtf_retro <- "/home/kamil/Desktop/n/snakemake_transposons2022/data/elegans_reference/ce11_rmsk_TE.nu.gtf"
gtf <- "/home/kamil/Desktop/n/snakemake_transposons2022/data/elegans_reference/c_elegans.PRJNA13758.WS283.canonical_geneset.gtf"

gtfretro <- fread(gtf_retro)
gtf1 <- fread(gtf)

gtf1 <- gtf1[gtf1$V3=="gene",]
gtf1 <- cbind(gtf1[,1], gtf1$V3,gtf1$V4,gtf1$V5,gtf1$V9)

midpoint_retro <- (gtfretro$V4+gtfretro$V5)/2
gtfretro <- cbind(gtfretro, midpoint_retro)
midpoint <- (gtf1$V3+gtf1$V4)/2
gtf1 <- cbind(gtf1, midpoint)

## distance analysis
full_matched <- list()
for(x in 1:dim(unique(gtf1[,1]))[1]){
  print(unique(gtf1[,1])[x])
  
  gtf1_chr1 <- gtf1[which(gtf1[,1]==paste(unique(gtf1[,1])[x])),]
  gtfretro_chr1 <- gtfretro[which(gtfretro[,1]==paste(unique(gtf1[,1])[x])),]
  if(paste(unique(gtf1[,1])[x])!="MtDNA"){
    nearest_gene <- vector()
    nearest_gene_midpoint <- vector()
    for(n in 1:length(gtfretro_chr1$midpoint_retro)){
      #print(n)
      rez <- abs(gtfretro_chr1$midpoint_retro[n]-gtf1_chr1$midpoint)
      dist <- min(rez)
      
      nearest_gene[n] <- gtf1_chr1$V5[which(rez==dist)[1]]
      nearest_gene_midpoint[n] <- gtf1_chr1$midpoint[which(rez==dist)[1]]
    }
    gtfretro_chr1 <- cbind(gtfretro_chr1,nearest_gene,nearest_gene_midpoint)
    full_matched[[x]] <- gtfretro_chr1
    }
  
}
result <- do.call(rbind,full_matched,)




