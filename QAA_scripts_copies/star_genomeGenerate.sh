#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=48GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=STAR_genomeGenerate_PS2    #optional: job name
#SBATCH --output=STAR_genomeGenerate_PS2_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=STAR_genomeGenerate_PS2_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

# be sure to change job, output, and error SBATCH names with every run.

# activate STAR environment
mamba activate bgmp_star

# run STAR genomeGenerate on Bi623 PS2:
/usr/bin/time -v STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index \
--genomeFastaFiles /projects/bgmp/joycew/bioinfo/Bi623/QAA/campylomormyrus.fasta \
--sjdbGTFfile /projects/bgmp/joycew/bioinfo/Bi623/QAA/campylomormyrus.gtf \
--genomeSAindexNbases 13

# run STAR genomeGenerate
# /usr/bin/time -v STAR \
# --runThreadN 8 \
# --runMode genomeGenerate \
# --genomeDir /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz/Danio_rerio.GRCz11.dna.ens114.STAR_2.7.11b \
# --genomeFastaFiles /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz/dre/Danio_rerio.GRCz11.dna.primary_assembly.fa \
# --sjdbGTFfile /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz/Danio_rerio.GRCz11.114.gtf \
