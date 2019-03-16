#!/usr/bin.bash

#########################################################################################
# DESCRIPTION
#
# Integrated assignment for CSHL CBW workshop (Day 3). Th assignment can be accessed at
# https://rnabio.org/module-08-appendix/0008/06/01/Integrated_Assignment/
#
#########################################################################################
#
# SCRIPT INFORMATION
#
# Author: Akshay Paropkari
# Date: 03/12/2019
#
#########################################################################################

# Create working directory to store all inputs and outputs

# We are already in $RNA_HOME, no need to `cd` into it. SKIP this code block.
#export RNA_HOME=/media/workspace/rnaseq  #~/workspace/rnaseq
#cd "$RNA_HOME" || echo -e "Could not change working directory to \"$RNA_HOME\""; exit
#echo "\n$(date "+%a %D %r"): We are in $(pwd)"

mkdir -p /media/workspace/rnaseq/integrated_assignment   #~/workspace/rnaseq/integrated_assignment/
export RNA_ASSIGNMENT=/media/workspace/rnaseq/integrated_assignment  #~/workspace/rnaseq/integrated_assignment/

# Set up environment variables
export RNA_DATA_DIR=$RNA_ASSIGNMENT/raw_reads
export RNA_REFS_DIR=$RNA_ASSIGNMENT/reference
export RNA_ILL_ADAPT=$RNA_ASSIGNMENT/adapter
export RNA_REF_INDEX=$RNA_REFS_DIR/Homo_sapiens.GRCh38
export RNA_REF_FASTA=$RNA_REF_INDEX.dna.primary_assembly.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf
export RNA_ALIGN_DIR=$RNA_ASSIGNMENT/hisat2

# Obtain reference, annotation, adapter and data files and place them in the integrated
# assignment directory
echo "$RNA_ASSIGNMENT" && cd "$RNA_ASSIGNMENT" || echo -e "Could not cd into \"$RNA_ASSIGNMENT\""
#ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/reference/
#ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/raw_reads/top_1mil/ raw_reads
#ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/adapter

echo "\n#########################################################################################"
echo "PART 0 : Obtaining Data and References"
echo "Goals:"
echo "Obtain the files necessary for data processing"
echo "Familiarize yourself with reference and annotation file format"
echo "Familiarize yourself with sequence FASTQ format"
echo "#########################################################################################"

echo "\n$(date "+%a %D %r")\nQ1.) How many items are there under the “reference” directory \(counting all files in all sub-directories\)? What if this reference file was not provided for you - how would you obtain/create a reference genome fasta file. How about the GTF transcripts file from Ensembl?"
echo "$(find reference/ -type f | wc -l)"

echo "\n$(date "+%a %D %r")\nQ2.) How many exons does the gene SOX4 have? How about the longest isoform of PCA3?"
grep "SOX4" $RNA_REF_GTF | awk '$3 == "exon"' | wc -l
grep "PCA3" $RNA_REF_GTF | awk '{print $5 - $4}' | sort -gr | head -1

echo "\n$(date "+%a %D %r")\nQ3.) How many samples do you see under the data directory?"
ls raw_reads/ | wc -l

echo "\n#########################################################################################"
echo "Part 1 : Data preprocessing"
echo "Goals:"
echo "Run quality check before and after cleaning up your data"
echo "Familiarize yourself with the options for Fastqc to be able to redirect your output"
echo "Perform adapter trimming on your data"
echo "Familiarize yourself with the output metrics from adapter trimming"
echo "#########################################################################################"

echo "\n$(date "+%a %D %r")\nQ4.) What metrics, if any, have the samples failed? Are the errors related?"
#mkdir -p $RNA_DATA_DIR/raw_fastqc
# fastqc $RNA_DATA_DIR/*.fastq.gz -o $RNA_DATA_DIR/raw_fastqc > $RNA_DATA_DIR/raw_fastqc/fastqc.log
echo "All samples fail Per base sequence content metric"

echo "\n$(date "+%a %D %r")\nQ5.) What average percentage of reads remain after adapter trimming? Why do reads get tossed out?"
mkdir -p $RNA_DATA_DIR/trimmed
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155055_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155055_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155055
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155056_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155056_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155056
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155057_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155057_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155057
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155058_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155058_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155058
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155059_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155059_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155059
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155060_1.fastq.gz  --reads2 $RNA_DATA_DIR/SRR7155060_2.fastq.gz --target $RNA_DATA_DIR/trimmed/SRR7155060

echo "\n$(date "+%a %D %r")\nQ6.) What sample has the largest number of reads after trimming?"

exit
echo "\n#########################################################################################"
echo "PART 2: Data alignment"
echo "Goals:"
echo "Familiarize yourself with HISAT2 alignment options"
echo "Perform alignments"
echo "Obtain alignment summary"
echo "Convert your alignment into compressed bam format"
echo "#########################################################################################"


echo "\n$(date "+%a %D %r")\nQ7.) How would you obtain summary statistics for each aligned file?"


echo "\n$(date "+%a %D %r")\nQ8.) Approximatly how much space is saved by converting the sam to a bam format?"



# Try viewing genes such as TP53 to get a sense of how the data is aligned. To do this:
# Load up IGV
# Change the reference genome to “Human hg38” in the top-left category
# Click on File > Load from URL, and in the File URL enter: “http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/transfected.bam”. Repeat this step and enter “http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/control.bam” to load the other bam, where ## is your student number for the AWS instance.
# Right-click on the alignments track in the middle, and Group alignments by “Library”
# Jump to TP53 by typing it into the search bar above


echo "\n$(date "+%a %D %r")\nQ9.) What portion of the gene do the reads seem to be piling up on? What would be different if we were viewing whole-genome sequencing data?"

echo "\n$(date "+%a %D %r")\nQ10.) What are the lines connecting the reads trying to convey?"


echo "\n#########################################################################################"
echo "PART 3: Expression Estimation"
echo "Goals:"
echo "Familiarize yourself with Stringtie options"
echo "Run Stringtie to obtain expression values"
echo "Obtain expression values for the gene SOX4"
echo "Create an expression results directory, run Stringtie on all samples, and store the results in appropriately named subdirectories in this results dir"
echo "#########################################################################################"


echo "\n$(date "+%a %D %r")\nQ11.) How do you get the expression of the gene SOX4 across the transfect and control samples?"


echo "\n#########################################################################################"
echo "PART 4: Differential Expression Analysis"
echo "Goals:"
echo "Perform differential analysis between the transfected and control samples"
echo "Check if is differentially expressed"
echo "#########################################################################################"

echo "\n$(date "+%a %D %r")\nQ12.) Are there any significant differentially expressed genes? How many in total do you see? If we expected SOX4 to be differentially expressed, why don’t we see it in this case?"

echo "\n$(date "+%a %D %r")\nQ13.) What plots can you generate to help you visualize this gene expression profile?"

