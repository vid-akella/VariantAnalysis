#!/bin/bash
##############################################################################
##############               BASH PIPELINE FOR                  ##############
##############                VARIANT CALLING                   ##############
##############################################################################
##############               author: Vidya Akella               ##############
##############               vidya.akella@gmail.com             ##############
##############################################################################

##############################################################################
##############      STEP 1: check if necessary tools are        ##############
##############           installed and install if absent        ##############
##############################################################################

# BWA MEM
if ! command -v bwa mem;
then 
    echo "Installing bwa!"
    sudo apt-get -y install bwa
fi

# SAMTOOLS
if ! command -v samtools;
then 
    echo "Installing samtools!"
    sudo apt-get -y install samtools
fi

# BCFTOOLS
if ! command -v bcftools;
then 
    echo "Installing bcftools!"
    sudo apt-get -y install bcftools
fi

##############################################################################
##############      STEP 2: Check for reference inputs          ##############
##############      or use preset references                    ##############
##############################################################################

# sample name is the first parameter passed to the script
sample_name=$1

# this block allows passing bwa ref to script as 2nd parameter; if not will use preset path 
if [ -z "$2" ]
then
    echo "No bwa ref input, using preset ref"
    bwa_ref="../refs/hg19.fa"
else
    bwa_ref=$2
fi 

# this block allows passing samtools ref to script as 3rd parameter; if not will use preset path 
if [ -z "$3" ]
then
    echo "No sam ref input, using preset ref"
    sam_ref="../refs/hg19.fa"
else
    sam_ref=$3
fi

##############################################################################
##############      STEP 3: Map reads to ref with BWA           ##############
##############                                                  ##############
##############################################################################

# run alignment with bwa
bwa mem $bwa_ref ${sample_name}1.fastq.gz ${sample_name}2.fastq.gz > $sample_name.sam

##############################################################################
##############      STEP 4: Generate sorted alignment           ##############
##############                                                  ##############
##############################################################################

# sort alignment with samtools
samtools sort ${sample_name}.sam > ${sample_name}_sorted.sam

##############################################################################
##############      STEP 5: Variant calling.                    ##############
##############                                                  ##############
##############################################################################
# variant calling with mpileup and bcftools
samtools mpileup -Ou -f ${sam_ref} ${sample_name}_sorted.sam | bcftools call -vmO v -o ${sample_name}.vcf


##############################################################################
##############      STEP 6: Variant call summaries              ##############
##############                                                  ##############
##############################################################################

#filter vcf for variant and mapping metrics and output in tabular format
grep -v '^#' NA12878-20k_R.vcf | awk '{n=split($8,a,/[;=]/); print $1, $2, $6, a[1], a[2], a[n-1], a[n]}'>variant_qual.txt

#calculate and print average variant quality
awk '{sum+=$3} END {print "AVERAGE VARIANT QUALITY SCORE = ", sum/NR}' variant_qual.txt 

# filter vcf table with thresholds for variant qual score, depth of coverage and mapping quality
cat variant_qual.txt| awk '$3>=30 && $5>=10 && $7>=30 {print $1, $2, $3, $4, $5, $6, $7}' >varQualFilt.txt
