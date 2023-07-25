library("stringr")

file_folder <- snakemake@input$in1
file_folder2 <- unlist(file_folder)
output_folder <- unlist(snakemake@output[[1]])
output_folder2 <- unlist(snakemake@output[[2]])

print(file_folder2)

## nyms
s1 <- str_split_fixed(file_folder2, "__", 2)[,2]
s2 <- str_split_fixed(s1, "/", 2)[,1]
print(s2)

names <- s2
filter_cutoff <- as.numeric(snakemake@params[[1]])
filter_cutoff2 <- as.numeric(snakemake@params[[2]])

print("entering loop1")
fl_list <- list()
for(n in 1:length(file_folder2)){
  print(n)
  c <- read.csv(paste0(file_folder2[n]), sep="\t")
  c <- cbind(rep(names[n], dim(c)[1]),c)
  names(c)[1] <- "sample"
  fl_list[[n]] <- c
}

fl_list2 <- do.call(rbind,fl_list)

## filter
active <- fl_list2$only+fl_list2$X3prunon
inactive <- fl_list2$passive_sense+fl_list2$passive_antisense+fl_list2$antisense
fl_filtered <- fl_list2[(active/inactive)>filter_cutoff, ]
fl_filtered <- fl_filtered[fl_filtered$only>filter_cutoff2, ]


## write
head(fl_list2)
head(fl_filtered)
write.csv(fl_list2, output_folder)
write.csv(fl_filtered, output_folder2)