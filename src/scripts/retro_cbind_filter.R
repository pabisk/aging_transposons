library("stringr")
file_folder <- snakemake@input$in1
genex_file_folder <- snakemake@input$in2
output_folder <- unlist(snakemake@output[[1]])
genex_output_folder <- unlist(snakemake@output[[2]])
read_nr <- as.numeric(snakemake@params[[1]])

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
# summarize the retro data

print("entering loop1")
fl_list <- list()
for(n in 1:length(file_folder2)){
  print(n)
  #print(paste0(file_folder2[n]))
  if(n==1){
    fl_list[[n]] <- read.csv(paste0(file_folder2[n]))
    fl_list[[n]] <- data.frame(fl_list[[n]][,1:3])
  }else{
    fl_list[[n]] <- read.csv(paste0(file_folder2[n]))
    nym <- names(fl_list[[n]])[3]
    fl_list[[n]] <- data.frame(fl_list[[n]][,3])
    names(fl_list[[n]]) <- nym
    
  }
}

# summarize the genex data
print("entering loop2")
genex_fl_list <- list()
for(n in 1:length(genex_file_folder2)){
  print(n)
  if(n==1){
    genex_fl_list[[n]] <- read.delim(paste0(genex_file_folder2[n]))
    #print(head(genex_fl_list[[n]], n=3))
    genex_fl_list[[n]] <- data.frame(genex_fl_list[[n]][,1:2])
   # print(head(genex_fl_list[[n]], n=3))
  }else{
    genex_fl_list[[n]] <- read.delim(paste0(genex_file_folder2[n]))
    nym <- names(genex_fl_list[[n]])[2]
    genex_fl_list[[n]] <- data.frame(genex_fl_list[[n]][,2])
    names(genex_fl_list[[n]]) <- nym
    
  }
}

print("genex_fl_list")
#print(head(genex_fl_list, n=3))
genex_fl_list2 <- do.call(cbind,genex_fl_list)

print("do.call(cbind,fl_list)")
fl_list2 <- do.call(cbind,fl_list)
fl_list3 <- fl_list2[rowSums(as.matrix(fl_list2[,3:dim(fl_list2)[2]]))>read_nr,] ## filter retro reads lower than a certain threshold

####
print("writing to file")
write.csv(fl_list3, output_folder)
write.table(genex_fl_list2, genex_output_folder, col.names=F, sep=",")