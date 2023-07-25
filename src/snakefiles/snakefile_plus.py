###############################
## given a data/SRR_Acc_List.txt download these SRRs to data/reads_fq
## SRR_Acc_List.txt is the file you can get from the SRA run selector
## on each line is one SRR ID

from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
import shutil
import yaml
#config = yaml.full_load("config.yml")
#configfile: "config.yaml" 

with open("config.yml", "r") as yamlfile:
    data = yaml.load(yamlfile, Loader=yaml.FullLoader)
    print("Read successful")

print("executing snakefile_plus.py")
print("config[aria3c_parallel_downloads] :")
print(config["aria3c_parallel_downloads"])
dlnr=config["aria3c_parallel_downloads"]
print(f'--dl_count={dlnr}')

mypath="data/reads_fq"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
df = pd.read_csv ('data/SRR_Acc_List.txt', header=None) ## otherwise it will consider the first entry to be a header

if(len(df)>0) and (len(onlyfiles) == 0):
    print("the following accessions were provided")
    print(df)
    print("running a perl script to download files of interest. Write yes to proceed, no to exit. Make sure no fastq (.fastq.gz) and .log files are in the snakemake main folder, as downloaded files will be temporarily stored there.")
    ask_user=input()
    if (ask_user == "yes"):
        print("running sra_download.pl")
        pipe = subprocess.call(["perl", "src/scripts/sra_download.pl" , "data/SRR_Acc_List.txt", f'--dl_count={dlnr}'], stdout=subprocess.PIPE)
        sourcepath='.'
        sourcefiles = os.listdir(sourcepath)
        destinationpath = 'data/reads_fq'
        for file in sourcefiles:
            if file.endswith('.fastq.gz'):
                shutil.move(os.path.join(sourcepath,file), os.path.join(destinationpath,file))
            if file.endswith('.log'):
                os.makedirs(os.path.dirname(os.path.join("results/logs",file)), exist_ok=True)
                shutil.move(os.path.join(sourcepath,file), os.path.join("results/logs",file))
else:
        print("data/SRR_Acc_List.txt empty and/or fastq files already present in data/reads_fq. Thus proceeding with default script.")

print(" ")
###############################
## clean up reads
## the underscore nomenclature is absolutely essential for
## PE: we have _1.filending and _2.filending
## SE: we have _1.filending
print("single end reads will be automatically cleaned up, i.e. _1.filending will be added")

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
mypath="data/reads_fq"

for file in onlyfiles:
    #print("file")
    #print(file)
    if file.endswith('_1.fastq.gz') or file.endswith('_1.fq.gz') or file.endswith('_2.fq.gz') or file.endswith('_2.fastq.gz'):
        print("no renaming and cleaning necessary")
    else:
        if file.endswith('.fastq.gz') or file.endswith('.fq.gz'):
            print("renaming your files, assuming they are unnamed SE fastq files, adding underscore _1")
            core_filename = file.split(".")[0]
            file_name_p2 = file.split(".")[1]
            file_name_p3 = file.split(".")[2]
            #print("core_filename")
            #print(core_filename)
            #print("file_name_p2")
            #print(file_name_p2)
            os.rename(mypath+"/"+file, mypath+"/"+core_filename+"_1"+"."+file_name_p2+"."+file_name_p3)
print(" ")
##############
## given a csv create a sequence dict
## file_name,group_name
## 2_R1_001,A
## only provide read 1
## Requirments: all files have to have unique names, dots only for file endings, double underscore only for groups
## in the samples.txt you use the core file name i.e. no file ending and no _1 or _2 suffix for the paired end reads

from os import listdir
from os.path import isfile, join
import pandas as pd
df = pd.read_csv ('data/samples.txt')
print("this is the df detailing the sample groups:")
print(df)

mypath="data/reads_fq"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

## only run this code if a df with sample names and groups actually exists
if(len(df)>0 and sum('__' in s for s in onlyfiles)==0):
    print("sample df is not empty and files have no prefix, proceeding to use groups from sample.txt")
    ## generate dictionary from file
    s=set(df["group"])
    a_dict = {}
    for groups_from_df in s:
        a_dict[groups_from_df] = set(df["file"][df["group"]==groups_from_df])


    ## now rename files according to the dictionary file
    counter=0
    for x in onlyfiles:
        counter=counter+1
        #print(counter)
        #print(x)
        n=x.count(df["ending"][0]) ## if we find a file with the proper ending continue
    
        if(n>0):
            core_filename = x.split(".")[0]
            size = len(core_filename)
            # Slice string to remove last 2 characters from string
            core_filename = core_filename[:size - 2]
    
            #print("testing")
            get_group=df["group"][df["file"]==core_filename].astype("string")
            get_group=[str(x) for x in get_group][0]
            os.rename(mypath+"/"+x,mypath+"/"+get_group+"__"+x)
    
    #print(a_dict)
    #print("end")
else:
    print("sample df is empty OR samples already have a group prefix, leave samples as is")
    print(" ")


##############
## given a list of files split at a given seperator and create a python dictionary from them
## if we are given a list of files WITHOUT a group separator, rename files, with file name = group name
## since the steps downstream of mutect2 collapse every analysis by group this should simply allow us to run every analysis on a per-file basis
## example files 1
## mouse55_2.fq.gz, mouse55_1.fq.gz -> rename -> mouse55__mouse55_2.fq.gz, mouse55__mouse55_1.fq.gz -> output vcf -> mouse55.vcf
## example files 2
## group_A__mouse55_2.fq.gz, group_A__mouse55_1.fq.gz -> NO rename -> -> output vcf -> group_A.vcf
from os import listdir
from os.path import isfile, join
mypath="data/reads_fq"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

a_dict = {}

counter=0
for x in onlyfiles:
    #print(counter)
    counter=counter+1
    
    ## no separator, no group case
    if len(x.split("__")) == 1:
        temp = x.split("__")[0]
        temp = temp.split(".")[0]
        
        size = len(temp)
        # Slice string to remove last 2 characters from string
        temp = temp[:size - 2]
        
        key, value = temp, temp
        #print("kv " + key, value)
        os.rename(mypath+"/"+x,mypath+"/"+temp+"__"+x)
        #print(mypath+"/"+x)
        #print(mypath+"/"+temp+"__"+x)
        
    ## has separator, has group case    
    if len(x.split("__")) == 2:    
        #print("has separator")
        key, value = x.split("__")
        #print(value.split("."))
        value = value.split(".")[0]
    
        size = len(value)
        # Slice string to remove last 2 characters from string
        value = value[:size - 2]
    
    if key in a_dict:
        if value in a_dict[key]:
            #print("redundant entry")
            n=0
        else:
            a_dict[key].append(value)
            #print("append")
    else:
        a_dict[key] = [value]
        #print("adding")

#print("exit")
#print(a_dict)  
#print(type(a_dict))  
GROUPS3=a_dict

print("this is the dictionary that defines your sample groups")
#print(GROUPS)
print(GROUPS3)
#print(SAMPLES_prep)
print(" ")

"""
## santize a list of int's
# eg cleanup 0.25 -> 0pp25 in our config so that we can use these strings to name our files
def cleanup_config(config_list):
    new_strings=list()
    for string in config_list:
        new_string = str(string)
        new_string = new_string.replace(".", "pp")
        new_strings.append(new_string)
    return(new_strings)
#print(cleanup_config(config["filter_f_score"]))

def reverse_cleanup_config(config_list):
    new_strings=list()
    for string in config_list:
        new_string = str(string)
        new_string = new_string.replace("pp", ".")
        new_strings.append(new_string)
    return(new_strings)
#print(cleanup_config(config["filter_f_score"]))
"""
