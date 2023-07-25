library("stringr")
#library("vroom")
#library("tibble")
library("data.table")
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
output_folder <- args[2]
reads_file <- args[3]
ele_file <- args[4]
ele_file_discarded <- args[5]
# use the cleaned data as input
clean_name <- str_split_fixed(str_split_fixed(basename(file),"\\.",2)[,1],"_",2)[,1]


print(file)
print(output_folder)
print(reads_file)

print("reading in file")
#x <- vroom::vroom(file)
x <- fread(file, header=F)
x_discard <- fread(ele_file_discarded, header=F)
#x <- x[28056:dim(x)[1],] #check whether this nr is correct

read_counts <- read.table(reads_file)
print("read_counts:")
print(as.numeric(read_counts$V1))

###################
## calculate nr of discarded transposons
discard <- fread(ele_file_discarded, header=F)
actual <- fread(ele_file, header=F)
discard_nr <- sum(discard$V2)
actual_nr <- sum(actual$V2)
discard_fraction <- discard_nr/(discard_nr+actual_nr)


###################

##if(str_count(x[1,1], "ENS") < 1){stop("ERROR")} ## ensmus ensg ensr rat etc
print("transposing")
transp <- x[-1,]
transp_discard <- x_discard[-1,]

transp2 <- data.frame(cbind(str_split_fixed(transp$V1,":",3),transp$V2))
transp2_discard <- data.frame(cbind(str_split_fixed(transp_discard$V1,":",3),transp_discard$V2))
#tsum <- transp2$X4

major_types <- data.frame(unique(transp2$X3), rep(0,length(unique(transp2$X3))))
for(n in 1:length(unique(transp2$X3))){
  major_types[n,2] <- sum(as.numeric(transp2[which(transp2$X3==unique(transp2$X3)[n]),]$X4))
}
major_types_discard <- data.frame(unique(transp2_discard$X3), rep(0,length(unique(transp2_discard$X3))))
for(n in 1:length(unique(transp2_discard$X3))){
  major_types_discard[n,2] <- sum(as.numeric(transp2_discard[which(transp2_discard$X3==unique(transp2_discard$X3)[n]),]$X4))
}


minor_types <- data.frame(unique(transp2$X2), rep(0,length(unique(transp2$X2))))
for(n in 1:length(unique(transp2$X2))){
  minor_types[n,2] <- sum(as.numeric(transp2[which(transp2$X2==unique(transp2$X2)[n]),]$X4))
}

minor_types_discard <- data.frame(unique(transp2_discard$X2), rep(0,length(unique(transp2_discard$X2))))
for(n in 1:length(unique(transp2_discard$X2))){
  minor_types_discard[n,2] <- sum(as.numeric(transp2_discard[which(transp2_discard$X2==unique(transp2_discard$X2)[n]),]$X4))
}

major_types_discard[,2] <- major_types_discard[,2]/major_types[,2]
minor_types_discard[,2] <-minor_types_discard[,2]/minor_types[,2]

names(major_types_discard) <- c("X1",clean_name);names(minor_types_discard) <- c("X1",clean_name);
names(major_types) <- c("X1",clean_name)
readdf <- data.frame("reads:",as.numeric(read_counts$V1)); names(readdf) <- c("X1",clean_name)
actual_nr <- data.frame("actual_nr:",as.numeric(actual_nr));  names(actual_nr) <- c("X1",clean_name)
discard_nr <- data.frame("discard_nr:",as.numeric(discard_nr)); names(discard_nr) <- c("X1",clean_name)
discard_fraction <- data.frame("discard_fraction:",as.numeric(discard_fraction)); names(discard_fraction) <- c("X1",clean_name)

print(dim(major_types_discard))  
print(dim(major_types)) 

print(dim(minor_types_discard))  
print(dim(minor_types)) 

names(minor_types) <- c("X1",clean_name)

print(major_types)
all_types <- rbind(major_types,minor_types)
all_types <- rbind(all_types,
                   readdf,
                   actual_nr,
                   discard_nr,
                   discard_fraction,
                   major_types_discard,
                   minor_types_discard)

print(tail(all_types))
write.csv(all_types,output_folder)
#write.csv(minor_types,"minor_types.csv")
#write.csv(major_types,"major_types.csv")