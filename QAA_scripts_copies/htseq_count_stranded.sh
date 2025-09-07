#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1. have more if writing multithreaded code
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu    #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=htseq_count_stranded            #optional: job name
#SBATCH --output=htseq_count_stranded_%j.out        #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=htseq_count_stranded_%j.err         #optional: file to store stderr from job, %j adds the assigned jobID

conda activate QAA

# run one at a time.

# /usr/bin/time -v htseq-count \
#     -f bam \
#     -r pos \
#     -s yes \
#     -t exon \
#     -i gene_id \
#     -m union \
#     ../STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam \
#     ../campylomormyrus.gtf \
#     > SRR25630308_stranded.count

/usr/bin/time -v htseq-count \
    -f bam \
    -r pos \
    -s yes \
    -t exon \
    -i gene_id \
    -m union \
    ../STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam \
    ../campylomormyrus.gtf \
    > SRR25630394_stranded.count