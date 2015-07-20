#!usr/bin/bash

# Abstract: This bash script concatenates all the R1 fastq gzip files and all
#           the R2 fastq gzip files. Also, creates and stores the md5 for each
#           concatenated R1/R2 file.

# Author: Akshay Paropkari

echo -e "\nSTART\nPlease enter output filename and press ENTER."
read -r name

FILES1=(*R1_001.fastq.gz)
for ((i=0; i<${#FILES1[@]}; i+=2 ))
do
    echo -e "\nR1\nConcatenating ${FILES1[i]} and ${FILES1[i+1]} to ${name}_R1.fastq.gz"
    cat "${FILES1[i]}" "${FILES1[i+1]}" > "$name"_R1.fastq.gz
    concat1="$name"_R1.fastq.gz
    echo -e "Generating MD5 for ${concat1}."
    md5 "$concat1" > md5_R1.txt
done

FILES2=(*R2_001.fastq.gz)
for ((i=0; i<${#FILES2[@]}; i+=2 ))
do
    echo -e "\nR2\nConcatenating ${FILES2[i]} and ${FILES2[i+1]} to ${name}_R2.fastq.gz"
    cat "${FILES2[i]}" "${FILES2[i+1]}" > "$name"_R2.fastq.gz
    concat2="$name"_R2.fastq.gz
    echo -e "Generating MD5 for ${concat2}."
    md5 "$concat2" > md5_R2.txt
done

echo -e "FINISH\n"
