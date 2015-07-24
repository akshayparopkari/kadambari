#!usr/bin/bash

# Abstract: Iterate through all sample files and get md5 sums for R1 and R2
#           files in the sample folder.

# Date: 07/24/2015

# Author: Akshay Paropkari

# Set up for loop in parent directory
for dir in */
do
    printf "%s\n" "$dir"    # print the folder name which is being processed
    # cd into each folder in the directory
    cd "$dir" || { echo "cd into data folder failed! Please check your working directory." ; exit 1 ; }

    printf 'START\n'

    FILE1=(*R1_001.fastq.gz)
    FILE2=(*R2_001.fastq.gz)

#   Creating md5 files
    printf "md5 for %s" "${FILE1[*]}\n"
    md5 "$FILE1" > R1_md5.txt
    printf "Done\n"

    printf "md5 for %s" "${FILE2[*]}\n"
    md5 "$FILE2" > R2_md5.txt
    printf "Done\n"

    printf "FINISH\n\n"

    # cd up into the parent directory
    cd ../ || { echo "cd up into the parent folder failed! Please check your working directory." ; exit 1 ; }

done

