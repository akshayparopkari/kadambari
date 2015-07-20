#!usr/bin/bash

# Abstract: This shell script iterates through all data folders and
#           uses cutadapt package for trimming paired end reads.

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

    printf "Trimming V1-V3 spacers in Plate 1 for %s and %s\n" "${FILE1[*]}" "${FILE2[*]}"
#     cutadapt -g 'ADAPTER_FWD' -G 'ADAPTER_REV' -o P1_R1_V1-V3_out.fastq -p P1_R2_V1-V3_out.fastq "${FILE1[*]}" "${FILE2[*]}" > plate1_V1-V3_summary.txt
    printf "Done!\n"

    printf "Trimming V4-V5 spacers in Plate 1 for %s and %s\n" "${FILE1[*]}" "${FILE2[*]}"
#     cutadapt -g 'ADAPTER_FWD' -G 'ADAPTER_REV' -o P1_R1_V4-V5_out.fastq -p P1_R2_V4-V5_out.fastq "${FILE1[*]}" "${FILE2[*]}" > plate1_V4-V5_summary.txt
    printf "Done!\n"

    printf "Trimming V1-V3 spacers in Plate 2 for %s and %s\n" "${FILE1[*]}" "${FILE2[*]}"
#     cutadapt -g 'ADAPTER_FWD' -G 'ADAPTER_REV' -o P2_R1_V1-V3_out.fastq -p P2_R2_V1-V3_out.fastq "${FILE1[*]}" "${FILE2[*]}" > plate2_V1-V3_summary.txt
    printf "Done!\n"

    printf "Trimming V4-V5 spacers in Plate 2 for %s and %s\n" "${FILE1[*]}" "${FILE2[*]}"
#     cutadapt -g 'ADAPTER_FWD' -G 'ADAPTER_REV' -o P2_R1_V4-V5_out.fastq -p P2_R2_V4-V5_out.fastq "${FILE1[*]}" "${FILE2[*]}" > plate2_V4-V5_summary.txt

    printf "FINISH\n\n"

    # cd up into the parent directory
    cd ../ || { echo "cd into parent folder failed! Please check your working directory." ; exit 1 ; }

done
