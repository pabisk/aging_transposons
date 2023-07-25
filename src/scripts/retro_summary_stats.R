library("stringr")
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
output_folder <- args[2]
reads_file <- args[3]
genex_files <- args[4]
clean_name <- str_split_fixed(str_split_fixed(basename(file),"\\.",2)[,1],"__",2)[,1]


print(file)
print(output_folder)
print(reads_file)

print("reading in file")
x <- read.table(file)

read_counts <- read.table(reads_file)
print("read_counts:")
print(as.numeric(read_counts$V1))

print("transposing")
transp <- x[-1,]

transp2 <- data.frame(cbind(str_split_fixed(transp$V1,":",3),transp$V2))


major_types <- data.frame(unique(transp2$X3), rep(0,length(unique(transp2$X3))))
for(n in 1:length(unique(transp2$X3))){
  major_types[n,2] <- sum(as.numeric(transp2[which(transp2$X3==unique(transp2$X3)[n]),]$X4))
}

minor_types <- data.frame(unique(transp2$X2), rep(0,length(unique(transp2$X2))))
for(n in 1:length(unique(transp2$X2))){
  minor_types[n,2] <- sum(as.numeric(transp2[which(transp2$X2==unique(transp2$X2)[n]),]$X4))
}
readdf <- data.frame("reads:",as.numeric(read_counts$V1))
names(major_types) <- c("X1",clean_name)
names(minor_types) <- c("X1",clean_name)
names(readdf) <- c("X1",clean_name)
print(major_types)
all_types <- rbind(major_types,minor_types)
all_types <- rbind(all_types, readdf)
print(tail(all_types))

write.csv(all_types,output_folder)
