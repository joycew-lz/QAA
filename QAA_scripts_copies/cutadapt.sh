#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1. have more if writing multithreaded code
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu    #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=cutadapt             #optional: job name
#SBATCH --output=cutadapt_%j.out        #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=cutadapt_%j.err         #optional: file to store stderr from job, %j adds the assigned jobID

conda activate QAA

# /usr/bin/time -v cutadapt \
    # -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    # -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    # -o SRR25630308_1.trimmed.fastq.gz \
    # -p SRR25630308_2.trimmed.fastq.gz \
    # SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz \
    # > cutadaapt_report_SRR25630308.txt

/usr/bin/time -v cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o SRR25630394_1.trimmed.fastq.gz \
    -p SRR25630394_2.trimmed.fastq.gz \
    SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz \
    > cutadaapt_report_SRR25630394.txt
