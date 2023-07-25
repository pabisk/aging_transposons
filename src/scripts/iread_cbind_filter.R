library("stringr")
print("running iread_cbind_filter.R")

file_folder <- snakemake@input$in1
file_folder2 <- unlist(file_folder)
output_folder <- unlist(snakemake@output[[1]])

print(file_folder2)

## nyms
s1 <- str_split_fixed(file_folder2, "__", 2)[,2]
s2 <- str_split_fixed(s1, "/", 2)[,1]
print(s2)

names <- s2
print(length(output_folder))

print("entering loop1")
fl_list <- list()
for(n in 1:length(file_folder2)){
  print(n)
  if(n==1){
    c <- read.csv(paste0(file_folder2[n]), sep="\t")[,1:2]
    colnames(c) <- c("genes",names[n])
  }else{
    c <- data.frame(read.csv(paste0(file_folder2[n]), sep="\t")[,2])
    colnames(c) <- c(names[n])
  }
  fl_list[[n]] <- c
  print(dim(fl_list[[n]]))
  print(names[n])
}
fl_list2 <- do.call(cbind,fl_list)

## write
print("fl_list2")
head(fl_list2, n=3)
write.csv(fl_list2, output_folder, row.names=F)