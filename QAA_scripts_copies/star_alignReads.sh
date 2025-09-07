#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=48GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=STAR_alignReads_PS2    #optional: job name
#SBATCH --output=STAR_alignReads_PS2_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=STAR_alignReads_PS2_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

# be sure to change job, output, and error SBATCH names with every run.

# activate STAR environment
mamba activate bgmp_star

# run STAR genomeGenerate on Bi623 PS2 (do one at a time):
# /usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
# --outFilterMultimapNmax 3 \
# --outSAMunmapped Within KeepPairs \
# --alignIntronMax 1000000 \
# --alignMatesGapMax 1000000 \
# --readFilesCommand zcat \
# --readFilesIn /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_1.trimmed.paired.fastq.gz \
#                 /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_2.trimmed.paired.fastq.gz \
# --genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index \
# --outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/SRR25630308_

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_1.trimmed.paired.fastq.gz \
                /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_2.trimmed.paired.fastq.gz \
--genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index \
--outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/SRR25630394_

# # run STAR alignReads
# /usr/bin/time -v STAR --runThreadN 8 --runMode alignRseads \
# --outFilterMultimapNmax 3 \
# --outSAMunmapped Within KeepPairs \
# --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
# --readFilesCommand zcat \
# --readFilesIn /projects/bgmp/shared/Bi621/dre_WT_ovar12_R1.qtrim.fq.gz /projects/bgmp/shared/Bi621/dre_WT_ovar12_R2.qtrim.fq.gz \
# --genomeDir /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz/Danio_rerio.GRCz11.dna.ens114.STAR_2.7.11b \
# --outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz/dre/WT_ovar12_