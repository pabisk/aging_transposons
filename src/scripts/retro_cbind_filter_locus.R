library("stringr")
library("data.table")
file_folder <- snakemake@input$in1
genex_file_folder <- snakemake@input$in2
genex_instance_folder <- snakemake@input$in3

output_folder <- unlist(snakemake@output[[1]])
genex_output_folder <- unlist(snakemake@output[[2]])
genex_instance_out_folder <- unlist(snakemake@output[[3]])

read_nr <- as.numeric(snakemake@params[[1]])

setDTthreads(as.numeric(snakemake@params[[2]]))


print("snakemake@input")
print(snakemake@input)

print("file_folder")
print(head(file_folder,n=6))

print("genex_file_folder")
print(head(genex_file_folder,n=6))


file_folder2 <- unlist(file_folder)
print("file_folder2")
print(head(file_folder2,n=6))

genex_file_folder2 <- unlist(genex_file_folder)
print("genex_file_folder2")
print(genex_file_folder2)

genex_instance_folder2 <- unlist(genex_instance_folder)
print("genex_file_folder2")
print(genex_instance_folder2)


# summarize the genex data
print("entering loop1 - gen expression data")
genex_fl_list <- list()
spl <- str_split_fixed(genex_file_folder2,"/",6)[,6]; spl <- str_split_fixed(spl,"__",2)[,1]
for(n in 1:length(genex_file_folder2)){
  #print(n)
  if(n==1){
    genex_fl_list[[n]] <- fread(paste0(genex_file_folder2[n]))
    print(head(genex_fl_list[[n]]))
    genex_fl_list[[n]] <- data.frame(genex_fl_list[[n]][,1:2])
    names(genex_fl_list[[n]]) <- c("genes", spl[n])
  }else{
    genex_fl_list[[n]] <- fread(paste0(genex_file_folder2[n]))
    # genex_fl_list[[n]] <- genex_fl_list[[n]][rowSums(as.matrix(genex_fl_list[[n]][,2:dim(genex_fl_list[[n]])[2]]))>read_nr,]
    genex_fl_list[[n]] <- data.frame(genex_fl_list[[n]][,2])
    names(genex_fl_list[[n]]) <- spl[n]
  }

  
}
#print(head(genex_fl_list, n=3))
genex_fl_list2 <- do.call(cbind,genex_fl_list)
genex_fl_list2 <- genex_fl_list2[rowSums(as.matrix(genex_fl_list2[,2:dim(genex_fl_list2)[2]]))>read_nr,]
write.csv(genex_fl_list2, genex_output_folder)


# summarize the retro SUMMARY data
print("entering loop2 - transposon summary data")
fl_list <- list()
spl <- str_split_fixed(file_folder2,"/",5)[,5]
spl <- str_split_fixed(spl,"__",2)[,1]
for(n in 1:length(file_folder2)){
  #print(n)
  #print(paste0(file_folder2[n]))
  if(n==1){
    fl_list[[n]] <- read.csv(paste0(file_folder2[n]))
    fl_list[[n]] <- data.frame(fl_list[[n]][,1:4])
    names(fl_list[[n]]) <- c("X", "type", "subtype",spl[n])
  }else{
    fl_list[[n]] <- read.csv(paste0(file_folder2[n]))
    # nym <- names(fl_list[[n]])[3]
    # print(nym)
    fl_list[[n]] <- data.frame(fl_list[[n]][,4])
    names(fl_list[[n]]) <- spl[n]
    
  }
}
print("do.call(cbind,fl_list)")
fl_list2 <- do.call(cbind,fl_list)
head(fl_list2)
fl_list3 <- fl_list2[rowSums(as.matrix(fl_list2[,4:dim(fl_list2)[2]]))>read_nr,] ## filter retro reads lower than a certain threshold
write.csv(fl_list3, output_folder)


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
    #genex_instance_list[[n]] <- read.delim(paste0(genex_instance_folder2[n]))
    genex_instance_list[[n]] <- fread(paste0(genex_instance_folder2[n]))
    genex_instance_list[[n]] <- data.frame(genex_instance_list[[n]][,2])
    # genex_instance_list[[n]] <- genex_instance_list[[n]][rowSums(as.matrix(genex_instance_list[[n]][,1:dim(genex_instance_list[[n]])[2]]))>read_nr,]
    names(genex_instance_list[[n]]) <- spl[n]
  }
}
print("genex_instance_list")
#print(head(genex_fl_list, n=3))
genex_instance_list2 <- do.call(cbind,genex_instance_list)
genex_instance_list2 <- genex_instance_list2[rowSums(as.matrix(genex_instance_list2[,2:dim(genex_instance_list2)[2]]))>read_nr,] ## filter retro reads lower than a certain threshold
write.csv(genex_instance_list2, genex_instance_out_folder)



####
#print("writing to file")

