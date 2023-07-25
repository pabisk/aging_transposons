library("stringr")
library("data.table")
options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
file <- args[1]
output_folder <- args[2]
reads_file <- args[3]
clean_name <- str_split_fixed(str_split_fixed(basename(file),"\\.",2)[,1],"_",2)[,1]
setDTthreads(1)
print(paste(R.Version()[c("major", "minor")], collapse = "."))

print("A")
print(file)
print(output_folder)
print(reads_file)

print("reading in file")
#x <- vroom::vroom(file)
x <- fread(file, header=F)
#x_discard <- fread(ele_file_discarded, header=F)
#x <- x[28056:dim(x)[1],] #check whether this nr is correct
print("B")
read_counts <- read.table(reads_file)
print("read_counts:")
print(as.numeric(read_counts$V1))

###################
## calculate nr of discarded transposons
#actual <- fread(xx, header=F)
head(x)
actual_nr_transposons <- sum(x$V2)
print("C")

###################

##if(str_count(x[1,1], "ENS") < 1){stop("ERROR")} ## ensmus ensg ensr rat etc
print("transposing")
transp <- x[-1,]
#transp_discard <- x_discard[-1,]

transp2 <- data.frame(cbind(str_split_fixed(transp$V1,":",4),transp$V2))


major_types <- data.frame(unique(transp2$X4), rep("major",length(unique(transp2$X4))), rep(0,length(unique(transp2$X4))) )
for(n in 1:length(unique(transp2$X4))){

  major_types[n,3] <- sum(as.numeric(transp2[which(transp2$X4==unique(transp2$X4)[n]),]$X5))
}

minor_types <- data.frame(unique(transp2$X3), rep("minor",length(unique(transp2$X3))), rep(0,length(unique(transp2$X3))))
for(n in 1:length(unique(transp2$X3))){
  minor_types[n,3] <- sum(as.numeric(transp2[which(transp2$X3==unique(transp2$X3)[n]),]$X5))
}

micro_types <- data.frame(unique(transp2$X2), rep("minor2",length(unique(transp2$X2))), rep(0,length(unique(transp2$X2))))
for(n in 1:length(unique(transp2$X2))){
  micro_types[n,3] <- sum(as.numeric(transp2[which(transp2$X2==unique(transp2$X2)[n]),]$X5))
}

names(major_types) <- c("X1",clean_name,clean_name)
names(micro_types) <- c("X1",clean_name,clean_name)
readdf <- data.frame("reads:",as.numeric(read_counts$V1), as.numeric(read_counts$V1)); names(readdf) <- c("X1",clean_name,clean_name)
actual_nr_transposons <- data.frame("actual_nr_transposons:",as.numeric(actual_nr_transposons), as.numeric(actual_nr_transposons));  names(actual_nr_transposons) <- c("X1",clean_name,clean_name)
#discard_nr <- data.frame("discard_nr:",as.numeric(discard_nr)); names(discard_nr) <- c("X1",clean_name)
#discard_fraction <- data.frame("discard_fraction:",as.numeric(discard_fraction)); names(discard_fraction) <- c("X1",clean_name)

#print(dim(major_types_discard))  
print(dim(major_types)) 

#print(dim(minor_types_discard))  
print(dim(minor_types)) 

names(minor_types) <- c("X1",clean_name,clean_name)

print(major_types)
all_types <- rbind(major_types,minor_types,micro_types)
all_types <- rbind(all_types,
                   readdf,
                   actual_nr_transposons)

print(tail(all_types))
write.csv(all_types,output_folder)
#write.csv(minor_types,"minor_types.csv")
#write.csv(major_types,"major_types.csv")
