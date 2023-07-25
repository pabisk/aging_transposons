library("stringr")
library("data.table")


genex_instance_folder <- snakemake@input$in3


genex_instance_out_folder <- unlist(snakemake@output[[1]])
print("A")

read_nr <- as.numeric(snakemake@params[[1]])

setDTthreads(as.numeric(snakemake@params[[2]]))
# cbind all files in this folder
# remove those with fewer than N reads


print("snakemake@input")
print(snakemake@input)


genex_instance_folder2 <- unlist(genex_instance_folder)
print("genex_file_folder2")
print(genex_instance_folder2)

print("B")

# do NOT summarize the genex data

# summarize the instance data
print("entering loop3 - transposon locus instance data")
genex_instance_list <- list()
spl <- str_split_fixed(genex_instance_folder2,"/",6)[,6]; spl <- str_split_fixed(spl,"__",2)[,1]
print(spl)

for(n in 1:length(genex_instance_folder2)){
  cat(n)
  if(n==1){
    genex_instance_list[[n]] <- fread(paste0(genex_instance_folder2[n]))
    genex_instance_list[[n]] <- data.frame(genex_instance_list[[n]][,1:2])
    names(genex_instance_list[[n]]) <- c("genes", spl[n])
    
  }else{
    genex_instance_list[[n]] <- fread(paste0(genex_instance_folder2[n]))
    genex_instance_list[[n]] <- data.frame(genex_instance_list[[n]][,2])
     names(genex_instance_list[[n]]) <- spl[n]
  }
}
print("genex_instance_list")

genex_instance_list2 <- do.call(cbind,genex_instance_list)

print("C")
print(genex_instance_out_folder)
write.csv(genex_instance_list2, genex_instance_out_folder)




####
#print("writing to file")

