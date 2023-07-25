#!/bin/bash


input=$1 
output=$2

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $input  > $output
