#!/bin/bash


save_space_yes_no=$1 
target_to_delete=$2
target_to_delete2=$3 
end=$4

echo $save_space_yes_no
echo $target_to_delete

if [[ $save_space_yes_no == 'yes' ]]; then
    echo "DELETING input to save space (experimental!)"

    if [[ $end == 'paired' ]]; then
            echo rm -v $target_to_delete
            rm -v $target_to_delete
            rm -v $target_to_delete2
    else
        echo rm -v $target_to_delete
        rm -v $target_to_delete
    fi
    #head $target_to_delete > ${target_to_delete}_2
    #rm $target_to_delete
    #mv ${target_to_delete}_2 $target_to_delete
    #rm ${target_to_delete}_2
else
    echo "not deleting any input (default settings)"
fi 
