# Bi623 PS2: RNA-seq Quality Assessment Assignment (QAA)

# Objective:
- Process electric organ and/or skeletal muscle RNA-seq reads for a future differential gene expression analysis
- Learn how to use existing tools for quality assessment and read trimming
- Learn how to compare quality assessments to those created by my own software
- Learn how to align and count reads
- Learn how to summarize important information in a high-level report
- PS2 works in tandem with PS4 for this pipeline:
    - [PS2] Assess quality of reads --> Trim reads --> Align reads to genome/transcriptome --> Count reads --> [PS4] Normalize reads counts --> Identify differentially expressed genes

# Git FORK assignment onto Talapas; keeping track of files...:
- Git fork the assignment, then git clone my fork into Talapas.
- ```cd /projects/bgmp/joycew/bioinfo/Bi623/QAA```
- Lab notebook: on personal computer, push up to Bi621_Lab_Notebook online git repo. Where I take all notes and record all of the necessary information for the final report.
- Rmd file: just instructions, on the QAA folder in Talapas. I can also view this on the QAA online repo.

# Dataset:
- Using 2 RNA-seq files from two different electric fish studies (PRJNA1005245 and PRJNA1005244)
- For all steps below, process the two libraries separately.
- Find data assignments here: `/projects/bgmp/shared/Bi623/PS2/QAA_data_Assignments.txt`
    ```
    SRR25630308     Joyce
    SRR25630394     Joyce
    ```

# Part 0: setting up the dataset:
- Refer to Bi623 ICA1.
- Download from NCBI SRA: SRR25630308, SRR25630394
    - ```cd /projects/bgmp/joycew/bioinfo/Bi623/QAA```
    - Using bash commands:
    ```{bash}
    conda activate SRA_env
    # already have sra-tools installed into that environment we created in ICA1

    prefetch SRR25630308
    prefetch SRR25630394
    ```

- Dump into FASTQ files
    ```{bash}
    fasterq-dump SRR25630308
    spots read      : 31,824,352
    reads read      : 63,648,704
    reads written   : 63,648,704

    # output SRR25630308_1.fastq, SRR25630308_2.fastq
    ```

    ```{bash}
    fasterq-dump SRR25630394
    spots read      : 5,885,243
    reads read      : 11,770,486
    reads written   : 11,770,486

    # output SRR25630394_1.fastq, SRR25630394_2.fastq
    ```

- Zip those files:
    - download pigz (into QAA)
    ```{bash}
    conda install -c conda-forge pigz
    ```

    - I tried starting an interactive session:
    ```{bash}
    srun --account=bgmp --partition=bgmp --time=1:00:00 --pty bash
    ```

    ```{bash}
    conda activate QAA

    pigz SRR25630308_*

    pigz SRR25630394_*
    ```

    - But it still took so long... so I wrote an sbatch script, ```gzip.sh```: (future #SBATCH in the shell scripts won't be recorded in this lab notebook.)

    ```
    #!/bin/bash
    #SBATCH --account=bgmp                    #REQUIRED: which account to use
    #SBATCH --partition=bgmp                  #REQUIRED: which partition to use
    #SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1. have more if writing multithreaded code
    #SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
    #SBATCH --mail-user=joycew@uoregon.edu    #optional: if you'd like email
    #SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
    #SBATCH --time=2:00:00                    #optional: set a time limit
    #SBATCH --job-name=gzip_files             #optional: job name
    #SBATCH --output=gzip_files_%j.out        #optional: file to store stdout from job, %j adds the assigned jobID
    #SBATCH --error=gzip_files_%j.err         #optional: file to store stderr from job, %j adds the assigned jobID

    conda activate QAA

    pigz SRR25630308_*.fastq SRR25630394_*.fastq
    ```



- We are processing this data for use in a future assignment, so please keep these files organized.

- Rename the files to reflect Species_sample_tissue_age/size_sample#_readnumber.fastq.gz. (I did this near the end of the assignment, so nothing else is named that way. :) )
    - mv SRR25630308_1.fastq.gz Cco_com124_EO_6cm_2_1.fastq.gz
    - mv SRR25630308_2.fastq.gz Cco_com124_EO_6cm_2_2.fastq.gz
    - mv SRR25630394_1.fastq.gz Crh_rhy115_EO_adult_1_1.fastq.gz
    - mv SRR25630394_2.fastq.gz Crh_rhy115_EO_adult_1_2.fastq.gz

# Part 1: Read quality score distributions

1. Create a new conda environment called ```QAA``` and install ```FastQC```, ```cutadapt```, and ```Trimmomatic```. Note the versions installed.

    ```bash
    conda create -n QAA
    conda activate QAA
    
    conda install -c bioconda fastqc=0.12.1
    fastqc --version
    0.12.1

    conda install -c bioconda cutadapt
    cutadapt --version
    2.6

    conda install -c bioconda trimmomatic
    trimmomatic -version
    0.40
    ```
2. Using FastQC via the command line on Talapas, produce plots of the per-base quality score distributions for R1 and R2 reads. Also, produce plots of the per-base N content, and comment on whether or not they are consistent with the quality score plots. Remember to /usr/bin/time.

    ```
    /usr/bin/time -v fastqc SRR25630308_1.fastq.gz
    /usr/bin/time -v fastqc SRR25630308_2.fastq.gz
    /usr/bin/time -v fastqc SRR25630394_1.fastq.gz
    /usr/bin/time -v fastqc SRR25630394_2.fastq.gz
    ```

    ```
    Command being timed: "fastqc SRR25630308_1.fastq.gz"
    User time (seconds): 192.26
    System time (seconds): 8.36
    Percent of CPU this job got: 99%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 3:21.07
    Maximum resident set size (kbytes): 686524

    Command being timed: "fastqc SRR25630308_2.fastq.gz"
    User time (seconds): 196.48
    System time (seconds): 9.35
    Percent of CPU this job got: 101%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 3:23.20

    Command being timed: "fastqc SRR25630394_1.fastq.gz"
    User time (seconds): 39.90
    System time (seconds): 1.92
    Percent of CPU this job got: 105%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:39.60

    Command being timed: "fastqc SRR25630394_2.fastq.gz"
    User time (seconds): 39.82
    System time (seconds): 2.04
    Percent of CPU this job got: 106%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:39.24
    Maximum resident set size (kbytes): 664552
    ```
    
    - This step asks me to compare the per-base N content plots and the per-base quality score plots for R1 and R2 of each sample, and see if they are consistent. If there is a drop in the per-base quality score, there should be a corresponding increase in N content in those same regions. Otherwise, if the quality scores are consistently high, the N content should remain very low and unchanging. Across both samples, and for both Reads 1 and 2 for SRR25630308 and SRR25630394, the per-base N content plots and per-base quality score plots are pretty consistent. Since both reads for both samples have very high quality scores overall, all the per-base N content plots look similar-- with a horizontal line at essentially 0 N content across all bases. This means there are very few gaps in the reads, and that the slight drops in quality scores as the base pair region increases across all reads are not because of gaps.

3. Run the quality score plotting script from the Demultiplexing assignment in Bi622. Describe how the FastQC quality score distribution plots compare to your own. If different, propose an explanation. Also, does the runtime differ? Mem/CPU usage? If so, why?

    - Check fastq.gz files for the read lengths for each of the reads.

        ```
        cd /projects/bgmp/joycew/bioinfo/Bi623/QAA

        zcat SRR25630308_1.fastq.gz | head
        zcat SRR25630308_2.fastq.gz | head
        zcat SRR25630394_1.fastq.gz | head
        zcat SRR25630394_2.fastq.gz | head

        # read lengths = 150 for all 4 reads.
        ```
    - In the QAA envrionment, ```conda install matplotlib```

    - Go to demultiplexing folder: ```cd /projects/bgmp/joycew/bioinfo/Bi622/Demultiplex/Assignment-the-first```. Edit ```part1_sbatch.sh``` in VSCode (didn't bother changing SBATCH job name, output, or error names as I ran one at a time). Run one at a time: ```sbatch part1_sbatch.sh```

        ```
        Command being timed: "./part1_script.py -f /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_1.fastq.gz -l 150"
        User time (seconds): 813.24
        System time (seconds): 0.22
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 13:36.12
        Maximum resident set size (kbytes): 62824

        Command being timed: "./part1_script.py -f /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_2.fastq.gz -l 150"
        User time (seconds): 831.04
        System time (seconds): 0.63
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 14:00.05
        Maximum resident set size (kbytes): 65068

        Command being timed: "./part1_script.py -f /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_1.fastq.gz -l 150"
        User time (seconds): 151.96
        System time (seconds): 0.16
        Percent of CPU this job got: 96%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:37.55
        Maximum resident set size (kbytes): 64876

        Command being timed: "./part1_script.py -f /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_2.fastq.gz -l 150"
        User time (seconds): 162.60
        System time (seconds): 0.22
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:43.27
        Maximum resident set size (kbytes): 62800
        ```

    - Rename the 4 png files. Move png files to ```cd /projects/bgmp/joycew/bioinfo/Bi623/QAA``` folder. I deleted the error and output files from the demultiplexing folder, from running these ps2 jobs, after pasting the usr time info, to keep the folders neat.

    - To view the FastQC quality score distribution plots: 
    
        ```
        unzip SRR25630308_1_fastqc.zip
        unzip SRR25630308_2_fastqc.zip
        unzip SRR25630394_1_fastqc.zip
        unzip SRR25630394_2_fastqc.zip
        ```

    Enter the fastqc folders for each of these and look at the **per_base_quality** images for the FastQC quality score distribution plots.

    - The FastQC quality score distribution plots are very similar compared to my own. I tried to track peaks and dips in quality scores, and for all both reads in both RNA-Seq samples, the quality score trends were similar and there were no discrepancies between which base pair positions in the reads experienced increases or decreases in quality scores. This makes sense, since both are using Phred quality scores from the same FASTQ files. The y-axis on my demultiplexing script is much more zoomed in for all plots, however, which allows the user to track read quality scores per base pair more efficiently, whereas the FASTQC quality score plots give the users an overarching idea of how good each read is by indicating where each base position sits in terms of high, mediocore, or low quality scores.

    - The runtimes between my demultiplexing script and FASTQC differ, as my script took longer on average for all four files. The FASTQC jobs took a higher percentage of CPU. My script, written in Python, was slower since it read through the FASTQ file line by line and stored scores in lists. The FASTQC command was faster but more computationally intensive since it keeps a lot more data structures in memory for additional quality report metrics beyond per-base qscore plots.
    
4. Comment on the overall data quality of your two libraries. Go beyond per-base qscore distributions. Examine the FastQC documentation for guidance on interpreting results and planning next steps. Make and justify a recommendation on whether these data are of high enough quality to use for further analysis.

- Looking at only the per-base qscore distributions, both libraries have high data quality. Going beyond these plots, I also looked at per-base sequence quality plots, and per-base sequence content plots to assess the overall data quality of the two libaries, SRR25630308 and SRR25630394. 
- For both read 1 and read 2 of the two libarries, the per-base sequence quality plots indicate most base pairs have a mean Phred quality score of around 33-36. This points to good quality reads that are usable for next steps, such as trimming, genome alignment, counting reads, and identifying differentially expressed genes. To confirm this, I also observed the per-base sequence content plots, which compares the proportion of the four DNA bases for each base position. For both the forward and reverse reads of SRR25630308 and SRR25630394, there is a lot of deviation between the proportion of A, T, C, and G in the first ~10 bp. Then, the four lines representing the DNA bases start mostly running in parallel with each other as the base pair number increases. This is expected of some types of libraries. RNA-Seq libraries in particular will always produce biased sequence composition near the start of the read (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html), and this is due to the fact that RNA-Seq libraries were produced by priming using random hexamers. These intrinsic biases in the positions at the start of the reads cannot be corrected by cutadapt trimming, but usually doesn't affect downstream analysis, so it's still safe to use these reads from the libraries SRR25630308 and SRR25630394 for the next steps of the differential gene expression analysis pipeline.
- 
# Part 2: Adapter trimming comparison

1. In the QAA environment, install Cutadapt and Trimmomatic and check the installations:

    ```
    conda install cutadapt 
    cutadapt --version (should be 5.0)
    4.4

    conda install trimmomatic=0.39

    trimmomatic -version (should be 0.39)
    ```

2. Using Cutadapt, properly trim adapter sequences from your assigned files. Be sure to read how to use Cutadapt. Use default settings.

    Try to determine what the adapters are on your own:
    - I couldn't find the adapters on my own.

    Use your Unix skills to search for the adapter sequences in your datasets and confirm the expected sequence orientations. Report the commands you used, the reasoning behind them, and how you confirmed the adapter sequences.
    
    - I used zgrep -c first to grep for the forward and reverse adapter sequences given to us by Hope. These results confirmed the adapter sequences as well as the expected sequence orientations, since I got reads when I searched for the forward adapter in the forward reads, and when I searched for the reverse adapter in the reverse reads, for both RNA-seq samples.

    ```
    # SRR25630308_1: R1 adapter
    zgrep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630308_1.fastq.gz
    2485944

    # SRR25630308_2: R2 adapter
    zgrep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" SRR25630308_2.fastq.gz  
    876639

    # SRR25630394_1: R1 adapter
    zgrep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630394_1.fastq.gz
    71909

    # SRR25630394_2: R2 adapter
    zgrep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" SRR25630394_2.fastq.gz  
    20218
    ```

    Create ```cutadapt.sh```, then ```sbatch cutadapt.sh```

    ```
    # run one at a time

    conda activate QAA

   /usr/bin/time -v cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o SRR25630308_1.trimmed.fastq.gz \
    -p SRR25630308_2.trimmed.fastq.gz \
    SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz \
    > cutadaapt_report_SRR25630308.txt

    /usr/bin/time -v cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o SRR25630394_1.trimmed.fastq.gz \
    -p SRR25630394_2.trimmed.fastq.gz \
    SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz \
    > cutadaapt_report_SRR25630394.txt
    ```

    ```
    Command being timed: "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o SRR25630308_1.trimmed.fastq.gz -p SRR25630308_2.trimmed.fastq.gz SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz"
	User time (seconds): 1031.28
	System time (seconds): 1.16
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 17:18.66
    Maximum resident set size (kbytes): 49548

    Command being timed: "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o SRR25630394_1.trimmed.fastq.gz -p SRR25630394_2.trimmed.fastq.gz SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz"
	User time (seconds): 178.57
	System time (seconds): 0.24
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:01.16
    Maximum resident set size (kbytes): 42120
    ```

    What proportion of reads (both R1 and R2) were trimmed?
    - For the SRR25630308 RNA-Seq sample:
        - Read 1 with adapter:               8,258,276 (25.9%)
        - Read 2 with adapter:               8,257,237 (25.9%)
    - For the SRR25630394 RNA-Seq sample:
        - Read 1 with adapter:                 458,160 (7.8%)
        - Read 2 with adapter:                 496,990 (8.4%)

    - Then, I used zgrep -c to grep for the forward and reverse adapter sequences after running cutadapt. As expected, these files have 0 counts of the adapter sequences because they were trimmed away.

    ```
    # SRR25630308_1 trimmed: R1 adapter
    zgrep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630308_1.trimmed.fastq.gz
    0

    # SRR25630308_2 trimmed: R2 adapter
    zgrep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" SRR25630308_2.trimmed.fastq.gz  
    0

    # SRR25630394_1 trimmed: R1 adapter
    zgrep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630394_1.trimmed.fastq.gz
    0

    # SRR25630394_2 trimmed: R2 adapter
    zgrep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" SRR25630394_2.trimmed.fastq.gz  
    0
    ```

    - I used zgrep -c first to grep for the forward and reverse adapter sequences given to us by Hope and saw several counts of R1 and R2 adapters in R1 and R2 files for both samples, respectively. As stated above, these results confirmed the adapter sequences as well as the expected sequence orientations, since I got reads when I searched for the forward adapter in the forward reads, and when I searched for the reverse adapter in the reverse reads, for both RNA-seq samples. After using zgrep -c again for the two adapter sequences after running cutadapt, however, I saw that these files have 0 counts of the adapter sequences, as cutadapt trimmed those away.


3. Use Trimmomatic to quality trim your reads (both R1 and R2). Specify the following, in this order:
- HEADCROP: 8 bases
- LEADING: quality of 3
- TRAILING: quality of 3
- SLIDING WINDOW: window size of 5 and required quality of 15
- MINLENGTH: 35 bases

    Create ```trimming.sh```, then ```sbatch trimming.sh```. Be sure to output compressed files and clear out all intermediate files (aka, we don't need the unpaired R1 and R2 files).

    ```
    # run one at a time

    conda activate QAA

    /usr/bin/time -v trimmomatic PE \
        -threads 8 \
        SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz \
        SRR25630308_1.trimmed.paired.fastq.gz SRR25630308_1.trimmed.unpaired.fastq.gz \
        SRR25630308_2.trimmed.paired.fastq.gz SRR25630308_2.trimmed.unpaired.fastq.gz \
        HEADCROP:8 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:5:15 \
        MINLEN:35

    /usr/bin/time -v trimmomatic PE \
        -threads 8 \
        SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz \
        SRR25630394_1.trimmed.paired.fastq.gz SRR25630394_1.trimmed.unpaired.fastq.gz \
        SRR25630394_2.trimmed.paired.fastq.gz SRR25630394_2.trimmed.unpaired.fastq.gz \
        HEADCROP:8 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:5:15 \
        MINLEN:35
    ```

    ```
    TrimmomaticPE: Started with arguments:
    -threads 8 SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz SRR25630308_1.trimmed.paired.fastq.gz SRR25630308_1.trimmed.unpaired.fastq.gz SRR25630308_2.trimmed.paired.fastq.gz SRR25630308_2.trimmed.unpaired.fastq.gz HEADCROP:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
    Quality encoding detected as phred33
    Input Read Pairs: 31824352 Both Surviving: 31472264 (98.89%) Forward Only Surviving: 200089 (0.63%) Reverse Only Surviving: 137970 (0.43%) Dropped: 14029 (0.04%)
    TrimmomaticPE: Completed successfully
        Command being timed: "trimmomatic PE -threads 8 SRR25630308_1.fastq.gz SRR25630308_2.fastq.gz SRR25630308_1.trimmed.paired.fastq.gz SRR25630308_1.trimmed.unpaired.fastq.gz SRR25630308_2.trimmed.paired.fastq.gz SRR25630308_2.trimmed.unpaired.fastq.gz HEADCROP:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35"
        User time (seconds): 1904.71
        System time (seconds): 35.69
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 32:30.08
        Maximum resident set size (kbytes): 914720
    
    TrimmomaticPE: Started with arguments:
    -threads 8 SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz SRR25630394_1.trimmed.paired.fastq.gz SRR25630394_1.trimmed.unpaired.fastq.gz SRR25630394_2.trimmed.paired.fastq.gz SRR25630394_2.trimmed.unpaired.fastq.gz HEADCROP:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
    Quality encoding detected as phred33
    Input Read Pairs: 5885243 Both Surviving: 5838224 (99.20%) Forward Only Surviving: 28057 (0.48%) Reverse Only Surviving: 17167 (0.29%) Dropped: 1795 (0.03%)
    TrimmomaticPE: Completed successfully
        Command being timed: "trimmomatic PE -threads 8 SRR25630394_1.fastq.gz SRR25630394_2.fastq.gz SRR25630394_1.trimmed.paired.fastq.gz SRR25630394_1.trimmed.unpaired.fastq.gz SRR25630394_2.trimmed.paired.fastq.gz SRR25630394_2.trimmed.unpaired.fastq.gz HEADCROP:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35"
        User time (seconds): 350.87
        System time (seconds): 6.14
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 5:59.96
        Maximum resident set size (kbytes): 987000
    ```


4. Plot the trimmed read length distributions for both **paired** R1 and **paired** R2 reads (on the same plot - yes, you will have to use Python or R to plot this. See ICA4 from Bi621). You can produce 2 different plots for your 2 different RNA-seq samples. There are a number of ways you could possibly do this. One useful thing your plot should show, for example, is whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you expect R1s and R2s to be adapter-trimmed at different rates and why.

    ```bash
    # Extract the sequence, then Measure its length, and Compute the distribution of lengths

    # For RNA-Seq Sample SRR25630308:

    # Read 1:
    zcat SRR25630308_1.trimmed.paired.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | sort | uniq -c | sort -n > SRR25630308_1.trimmed.paired_dist

    # Read 2:
    zcat SRR25630308_2.trimmed.paired.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | sort | uniq -c | sort -n > SRR25630308_2.trimmed.paired_dist

    # For RNA-Seq Sample SRR25630394:

    # Read 1:
    zcat SRR25630394_1.trimmed.paired.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | sort | uniq -c | sort -n > SRR25630394_1.trimmed.paired_dist

    # Read 2:
    zcat SRR25630394_2.trimmed.paired.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | sort | uniq -c | sort -n > SRR25630394_2.trimmed.paired_dist
    ```

    Created R script on my personal computer```trimmed_read_lengths_dist.rmd``` to plot.
    - First, download these trimmed paired distribution files from Talapas onto my personal computer: ```cd /Users/joycewang/bioinfo/Bi623/PS/PS2```
    ```scp joycew@login4.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/*.trimmed.paired_dist .```
    - Then, loaded in data.
    - Then, plotted in R.


    Upload this R script and the two png bar charts for the two samples, ```SRR25630308_read_length_dist.png``` and ```SRR25630394_read_length_dist.png```, from personal computer to Talapas: ```scp *.Rmd joycew@login4.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/``` and ```scp *.png joycew@login4.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/```

    Discuss whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you expect R1s and R2s to be adapter-trimmed at different rates and why.
    - The read length histograms for each RNA-Seq sample show nearly overlapping R1 and R2 distributions—I avoided having the plots completely overlap by using position = "dodge" in ggplot2. This indicates that neither read direction was trimmed significantly more than the other. Both samples retained similar overall shapes with for both R1 and R2, suggesting minimal adapter contamination in either R1 or R2 for SRR25630308 and SRR25630394. This outcome is expected because the insert sizes were likely long enough relative to the sequencing read length (150 bp) for both R1 and R2, so that adapters appeared and were trimmed at similar rates in both read directions. After adapter trimming, the read length distributions reflect the biological sequences without bias from leftover adapter sequences.
    
5. Bonus - Run FastQC on your trimmed data. Comment on differences you observe between the trimmed and untrimmed data. Include any figures needed to support your conclusions.

# Part 3: Alignment to a genome and strand-specificity
- Align reads to the genome (not a transcriptome) next.
- It also gives information on quality:
    - Contamination: what percent of the reads map to the genome/transcriptome?
- Aligner: STAR: Spliced Transcripts Alignment to a Reference (align RNA-seq reads to a reference genome!)
- Evaluator of alignment: Picard

1. Install things
```
conda activate QAA

conda install star
conda install picard=2.18
conda install samtools
conda install numpy # all requested packages already installed
conda install matplotlib # all requested packages already installed
conda install htseq

```
2. Copy Campylomormyrus compressirostris genome fasta and gff file: ```cd /projects/bgmp/joycew/bioinfo/Bi623/QAA```, ```cp /projects/bgmp/shared/Bi623/PS2/campylomormyrus.fasta .```, ```cp /projects/bgmp/shared/Bi623/PS2/campylomormyrus.gff .```

3. Align the reads to your C. compressirostris database using a splice-aware aligner. Use the settings specified in PS8 from Bi621...

- convert GFF file into GTF file:
    ```
    conda install bioconda::gffread

    gffread campylomormyrus.gff -T -o campylomormyrus.gtf
    ```

- use the PS8 scripts now to 1) generate genome indexes and 2) map/align reads to the genome.

    ```
    cd /projects/bgmp/joycew/bioinfo/Bi621/PS/ps8-joycew-lz
    ```

- generate genome indexes:
    - first, ```mkdir STAR_index``` for the outputs to go into in the QAA directory.
    - in ```star_genomeGenerate.sh```, edit the SBATCH and comment out the previous commands, and input the following. Run with ```sbatch star_genomeGenerate.sh```

    ```
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
    ----genomeSAindexNbases 13 # added after warning
    ```

    ```
    !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=862592683, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13
    ```
    - Reran with suggestion.

    ```
    Command being timed: "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index --genomeFastaFiles /projects/bgmp/joycew/bioinfo/Bi623/QAA/campylomormyrus.fasta --sjdbGTFfile /projects/bgmp/joycew/bioinfo/Bi623/QAA/campylomormyrus.gtf --genomeSAindexNbases 13"
	User time (seconds): 1528.03
	System time (seconds): 21.59
	Percent of CPU this job got: 456%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 5:39.37
    Maximum resident set size (kbytes): 25219012
    ```

- align data:
    - first, ```mkdir STAR_align``` for the outputs to go into in the QAA directory.
    - in ```star_alignReads.sh```, edit the SBATCH and comment out the previous commands, and input the following. Run with ```sbatch star_alignReads.sh```

    ```
    #SBATCH --job-name=STAR_alignReads_PS2    #optional: job name
    #SBATCH --output=STAR_alignReads_PS2_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
    #SBATCH --error=STAR_alignReads_PS2_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

    # be sure to change job, output, and error SBATCH names with every run.

    # activate STAR environment
    mamba activate bgmp_star

    # run STAR genomeGenerate on Bi623 PS2 (do one at a time):
    /usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_1.trimmed.paired.fastq.gz \
                    /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_2.trimmed.paired.fastq.gz \
    --genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index \
    --outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/SRR25630308_

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
    ```

    ```
    Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_1.trimmed.paired.fastq.gz /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630308_2.trimmed.paired.fastq.gz --genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index --outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/SRR25630308_"
	User time (seconds): 4661.50
	System time (seconds): 14.35
	Percent of CPU this job got: 768%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 10:08.53
    Maximum resident set size (kbytes): 10236140

    Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_1.trimmed.paired.fastq.gz /projects/bgmp/joycew/bioinfo/Bi623/QAA/SRR25630394_2.trimmed.paired.fastq.gz --genomeDir /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_index --outFileNamePrefix /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/SRR25630394_"
	User time (seconds): 700.87
	System time (seconds): 5.60
	Percent of CPU this job got: 704%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:40.23
    Maximum resident set size (kbytes): 10069892
    ```

4. Remove PCR duplicates using Picard MarkDuplicates. You may need to sort your reads with samtools before running Picard. ```cd /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align```
- Use samtools to sort the reads:
    ```
    # convert SAM to BAM:

    samtools view -bS SRR25630308_Aligned.out.sam -o SRR25630308_Aligned.bam

    samtools view -bS SRR25630394_Aligned.out.sam -o SRR25630394_Aligned.bam

    # sort BAM:
    samtools sort SRR25630308_Aligned.bam -o SRR25630308_Aligned.sorted.bam

    samtools sort SRR25630394_Aligned.bam -o SRR25630394_Aligned.sorted.bam

    ```
- Create ```/projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_align/picard.sh``` in the QAA/STAR_align directory and deduplicate PCR dupicates with Picard. Run with ```sbatch picard.sh```.
```
#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1. have more if writing multithreaded code
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu    #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=picard_dedup             #optional: job name
#SBATCH --output=picard_dedup_%j.out        #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=picard_dedup_%j.err         #optional: file to store stderr from job, %j adds the assigned jobID

conda activate QAA

# run one at a time.

picard MarkDuplicates \
    INPUT=SRR25630308_Aligned.sorted.bam \
    OUTPUT=SRR25630308_Aligned.sorted.dedup.bam \
    METRICS_FILE=SRR25630308.metrics \
    REMOVE_DUPLICATES=TRUE \
    VALIDATION_STRINGENCY=LENIENT

picard MarkDuplicates \
    INPUT=SRR25630394_Aligned.sorted.bam \
    OUTPUT=SRR25630394_Aligned.sorted.dedup.bam \
    METRICS_FILE=SRR25630394.metrics \
    REMOVE_DUPLICATES=TRUE \
    VALIDATION_STRINGENCY=LENIENT
```

- I didn't record user/bin/time for the picard MarkDuplicates removing PCR duplicates step.

- After all this, I renamed the ```STAR_Align``` directory to ```STAR_Align_and_picard```, so I know that I did the picard PCR deduplicating step there.

5. Using your script from PS8 in Bi621, report the number of mapped and unmapped reads from each of your **2 sam files post deduplication** with picard. Make sure that your script is looking at the bitwise flag to determine if reads are primary or secondary mapping (update/fix your script if necessary).
- I only have 2 BAM files post-deduplication. I need to convert BAM to SAM...
```
# Convert BAM to SAM:
# -h to include the header lines

cd /projects/bgmp/joycew/bioinfo/Bi623/QAA/STAR_Align_and_picard

samtools view -h SRR25630308_Aligned.sorted.dedup.bam > SRR25630308_Aligned.sorted.dedup.sam

samtools view -h SRR25630394_Aligned.sorted.dedup.bam > SRR25630394_Aligned.sorted.dedup.sam

``` 

- Since the instructions indicate I might be adjusting the script, I actually made a copy of ```ps8_script.py``` into ```QAA``` and called it ```report_mapped_unmapped.py```.
- All I ended up changing was the ```with open("...", 'r') as f:```.
- Run using ```./report_mapped_unmapped.py```, after ```chmod 755 report_mapped_unmapped.py```. Run twice, replacing the file name for both RNA-Seq samples in the code... 

- ```"/STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam"```:
    - **Number of MAPPED reads: 21839114**
    - **Number of UNMAPPED reads: 19106563**

- ```"/STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam"```: 
    - **Number of MAPPED reads: 6672478**
    - **Number of UNMAPPED reads: 1302451**


6. Count deduplicated reads that map to features using htseq-count. You should run htseq-count twice: once with --stranded=yes and again with --stranded=reverse. Use default parameters otherwise. You may need to use the -i parameter for this run.

    - Yes and reverse mean it’s a stranded library prep method...

    - ```mkdir htseq_count```.

    - make 2 sbatch scripts, one for stranded=yes, one for stranded=reverse. I'm technically running htseq 4 times (once for each sample, for both stranded=yes and stranded=reverse).

    - template: (also refer to documentation, https://htseq.readthedocs.io/en/release_0.11.1/count.html)
    ```
    htseq-count –s yes –-t exon -i gene_id -m union 
    path_to_samfile \ # NOTE: SAM file, not BAM
    path_to_gtf > outfile

    • S: Stranded parameter is based on library prep
    • Type: exon parameter restricts counting reads aligned to exons
    • Exon is the default, you can set this flag to remember this
    • i: gene_id sums counts across all features within a gene
    • Mode: Union mode restricts count of reads aligned to 1 gene (union is also the default, can set for your bookkeeping)
    ```

- stranded=yes: ```htseq_count_stranded.sh```, run with ```sbatch htseq_count_stranded.sh```

```
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

/usr/bin/time -v htseq-count \
    -f bam \
    -r pos \
    -s yes \
    -t exon \
    -i gene_id \
    -m union \
    ../STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam \
    ../campylomormyrus.gtf \
    > SRR25630308_stranded.count

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
```

```
Warning: Mate records missing for 34117 records; first such record: <SAM_Alignment object: Paired-end read 'SRR25630308.193008', not aligned>.
20489897 alignment record pairs processed.
	Command being timed: "htseq-count -f bam -r pos -s yes -t exon -i gene_id -m union ../STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam ../campylomormyrus.gtf"
	User time (seconds): 1225.98
	System time (seconds): 9.82
	Percent of CPU this job got: 99%

Warning: Mate records missing for 3437 records; first such record: <SAM_Alignment object: Paired-end read 'SRR25630394.233', not aligned>.
3989183 alignment record pairs processed.
	Command being timed: "htseq-count -f bam -r pos -s yes -t exon -i gene_id -m union ../STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam ../campylomormyrus.gtf"
	User time (seconds): 274.43
	System time (seconds): 1.00
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:39.94

```

- stranded=reverse: ```htseq_count_reverse.sh```, run with ```sbatch htseq_count_reverse.sh```

```
#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1. have more if writing multithreaded code
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu    #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --time=2:00:00                    #optional: set a time limit
#SBATCH --job-name=htseq_count_reverse            #optional: job name
#SBATCH --output=htseq_count_reverse_%j.out        #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=htseq_count_reverse_%j.err         #optional: file to store stderr from job, %j adds the assigned jobID

conda activate QAA

# run one at a time.

/usr/bin/time -v htseq-count \
    -f bam \
    -r pos \
    -s reverse \
    -t exon \
    -i gene_id \
    -m union \
    ../STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam \
    ../campylomormyrus.gtf \
    > SRR25630308_reverse.count

/usr/bin/time -v htseq-count \
    -f bam \
    -r pos \
    -s reverse \
    -t exon \
    -i gene_id \
    -m union \
    ../STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam \
    ../campylomormyrus.gtf \
    > SRR25630394_reverse.count
```

```
Warning: Mate records missing for 34117 records; first such record: <SAM_Alignment object: Paired-end read 'SRR25630308.193008', not aligned>.
20489897 alignment record pairs processed.
	Command being timed: "htseq-count -f bam -r pos -s reverse -t exon -i gene_id -m union ../STAR_Align_and_picard/SRR25630308_Aligned.sorted.dedup.sam ../campylomormyrus.gtf"
	User time (seconds): 1226.20
	System time (seconds): 4.20
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 20:34.13
    Maximum resident set size (kbytes): 229508

Warning: Mate records missing for 3437 records; first such record: <SAM_Alignment object: Paired-end read 'SRR25630394.233', not aligned>.
3989183 alignment record pairs processed.
	Command being timed: "htseq-count -f bam -r pos -s reverse -t exon -i gene_id -m union ../STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam ../campylomormyrus.gtf"
	User time (seconds): 275.65
	System time (seconds): 1.09
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:41.47
    Maximum resident set size (kbytes): 181936
```
- the outputs are:

```
cd /projects/bgmp/joycew/bioinfo/Bi623/QAA/htseq_count

SRR25630308_reverse.count
SRR25630308_stranded.count
SRR25630394_reverse.count
SRR25630394_stranded.count
```
7. Demonstrate convincingly whether or not the data are from "strand-specific" RNA-Seq libraries and which stranded= parameter should you use for counting your reads for a future differential gene expression analyses. Include any commands/scripts used. Briefly describe your evidence, using quantitative statements (e.g. "I propose that these data are/are not strand-specific, because X% of the reads are y, as opposed to z."). 

- The kit used during library preparation for the RNA-Seq libraries of both SRR25630308 and SRR25630394 is called "NEXTFLEX Rapid Directional RNA-seq 2.0", and it is considered a more modern technique for RNA-seq. The company website states that the libraries produced with this technology are stranded, and that dUTP is incorporated during the second strand synthesis step. Since the strand containing dUTP is not amplified during library amplification, this kit also retains the original strands' orientations (Source: https://www.revvity.com/product/nex-rapid-dir-rna-seq-kit-2-0-8rxn-nova-5198-01). This determines that the data we are using are strand-specific RNA-Seq libraries. Also, based on the paper provided in the instructions by Krishna et al., kits that use first strand cDNA libraries that mark dUTP in the second strand require  ```--s reverse``` for the ```htseq-count``` tool's stranded option (Source: Krishna A Srinivasan, Suman K Virdee, Andrew G McArthur, Strandedness during cDNA synthesis, the stranded parameter in htseq-count and analysis of RNA-Seq data, Briefings in Functional Genomics, Volume 19, Issue 5-6, September-November 2020, Pages 339–342, https://doi.org/10.1093/bfgp/elaa010). Thus, we should use the ```stranded=reverse``` parameter for counting the reads for future differential gene expression analyses.
- To prove this, we are going to count the reads from all four files (each sample with htseq-count ran with stranded=reverse AND stranded=yes) and compare the percentage of reads that mapped to a feature (gene) for each file, using code from Bi621's ICA4.

    ```
    SRR25630308_reverse.count
    SRR25630308_stranded.count
    SRR25630394_reverse.count
    SRR25630394_stranded.count
    ```
## (refer to Bi621 ICA4) Determine the percentage of reads that mapped to a feature (gene) for each file, separately.

1. Sum the number of reads that mapped to a feature. (one line command for each file):
    ```
    cat SRR25630308_reverse.count | grep -v "^__" | awk '{sum += $2} END {print sum}'
    cat SRR25630308_stranded.count | grep -v "^__" | awk '{sum+=$2} END {print sum}'
    ```
    ```
    Number of reads mapped, SRR25630308_reverse: 6308181
    Number of reads mapped, SRR25630308_stranded: 347629
    ```

    ```
    cat SRR25630394_reverse.count | grep -v "^__" | awk '{sum += $2} END {print sum}'
    cat SRR25630394_stranded.count | grep -v "^__" | awk '{sum+=$2} END {print sum}'
    ```
    ```
    Number of reads mapped, SRR25630394_reverse: 2068866
    Number of reads mapped, SRR25630394_stranded: 115361
    ```
2. Calculate the total number of reads. (one line command for each file)

    ```
    awk '{sum += $2} END {print sum}' SRR25630308_reverse.count
    awk '{sum += $2} END {print sum}' SRR25630308_stranded.count
    ```

    ```
    Total reads, SRR25630308_reverse: 20489897
    Total reads, SRR25630308_stranded: 20489897
    ```

    ```
    awk '{sum += $2} END {print sum}' SRR25630394_reverse.count
    awk '{sum += $2} END {print sum}' SRR25630394_stranded.count
    ```

    ```
    Total reads, SRR25630394_reverse: 3989183
    Total reads, SRR25630394_stranded: 3989183
    ```

3. Determine the percentage of reads mapped by dividing the number of mapping reads by the total number of reads. (feel free to use a calculator, remember to report a percentage)

    ```
    For SRR25630308:
    Percent of reads mapped, SRR25630308_reverse: 6308181/20489897 x (100%) = 30.8%
    Percent of reads mapped, SRR25630308_stranded: 347629/20489897 x (100%) = 1.8%
    ```

    ```
    For SRR25630394:
    Percent of reads mapped, SRR25630394_reverse: 2068866/3989183 x (100%) = 51.9%
    Percent of reads mapped, SRR25630394_stranded: 115361/3989183 x (100%) = 2.9%
    ```


4. Which file had more reads mapping to features? The numbers you're computing in the questions above may differ between the stranded=yes and stranded=reverse files.
    - For SRR25630308, stranded=reverse assigned 30.8% of reads to features while stranded=yes assigned 1.8%. For SRR25630394, stranded=reverse assigned 51.9% vs. stranded=yes 2.9%. That means for both SRR25630308 and SRR25630394 RNA-seq libraries, the count files that had more reads mapping to features (or genes) was where htseq-count used the stranded=reverse option.
    - I propose that these data are strand-specific, because the majority of reads mapped to features only when using stranded=reverse. This large discrepancy indicates that the libraries are directional and affirms my earlier decision that the read counts file we should use are the ones where the stranded parameter is set to reverse.
    - We can then use SRR25630308_reverse.count and SRR25630394_reverse.count PS4's downstream analysis, where we will normalize the reads counts and identify differentially expressed genes in a DEG analysis.

- Copy the counts files that would be used in a future differential RNA-seq analysis from the ```htseq_count``` directory into a new ```htseq_count_files``` directory that will be submitted! Renamed to follow  the naming conventions from the RNAseq meta data file: YourSubmittedHTSeqCountsFileName column.
    - Made copies of ```SRR25630308_reverse.count``` and  ```SRR25630394_reverse.count``` into the ```htseq_count_files``` directory.
    - ```SRR25630308_reverse.count``` --> ```Cco_com124_EO_6cm_2_htseqcounts_[forORrev]stranded.txt```.
    - ```SRR25630394_reverse.count``` --> ```Crh_rhy115_EO_adult_1_htseqcounts_[forORrev]stranded.txt```



8. BONUS - Turn your commands from part 1 and 2 into a script with a loop going through your two SRA files

9. Remember, earlier, we were going to rename files?
- Rename the files to reflect Species_sample_tissue_age/size_sample#_readnumber.fastq.gz. (I did this near the end of the assignment, so nothing else is named that way. :) )
    ``` 
    mv SRR25630308_1.fastq.gz Cco_com124_EO_6cm_2_1.fastq.gz
    mv SRR25630308_2.fastq.gz Cco_com124_EO_6cm_2_2.fastq.gz
    mv SRR25630394_1.fastq.gz Crh_rhy115_EO_adult_1_1.fastq.gz
    mv SRR25630394_2.fastq.gz Crh_rhy115_EO_adult_1_2.fastq.gz
    ```
10. Go back and answer any unanswered questions. After this is done... create a PDF report.

# PDF report details:
- use Rmarkdown to create a report, ```QAA_report.pdf``` on my personal computer. Include:
    - all requested plots
    - answers to questions
    - mapped/unmapped read counts from PS8 script, in a nicely formatted table
    - submit this to the top level of this PS2 repo!
    - submit this to Canvas!

- Put plot photos used for this report into a separate folder, ```QAA_report_plots```
- Upload ```QAA_report.pdf``` and the images folder ```QAA_report_plots``` from my personal computer onto Talapas, then onto Git:
    ```
    cd ~/bioinfo/Bi623/PS/PS2

    scp QAA_report* joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/

    # to copy the images directory and all its content:
    scp -r QAA_report_plots joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/
    ```

# Talapas batch script/code organization --> QAA_scripts_copies:
- Manually MAKE A COPY OF EVERY script I used, and shove it in the ```/QAA/QAA_scripts_copies``` folder. Even if I made the script just for this assignment, or even if I made a copy from a script from another assignment into my ```QAA``` folder, or a folder nestled in the QAA folder, just COPY EVERYTHING into ```QAA_scripts_copies```.
- I'm doing this because, even though I used absolute paths instead of relative paths for most of my scripts, I still want to be careful and don't want to mess things up according to how I did/recorded things in my lab notebook... it's too much of a hassle.
- **THERE IS NO GUARANTEE THESE COPIES OF PY/BATCH SCRIPTS WILL WORK because of the paths used...**
- **Notes for next time:** just make a copy of any old scripts I need, and shove **all** scripts into a specific script folder (instead of out in the parent directory, or nestled in another subdirectory).
- I double checked that I made a copy of every .py and .sh script by cmd+F for ".py" and ".sh" in this lab notebook, and making sure I copied over every single one of these files mentioned in this lab notebook into ```QAA_scripts_copies```.

# Upload:
- FastQC plots: CHECK
    - in a folder with plots used for QAA report, ```QAA_report_plots```
    - quality is kind of bad in the way I organized it into one file for each read using Google Slides (for the convenience of putting together ```QAA_report.rmd``` and ```QAA_report.pdf```)...
    - so I made a copy of every plot's svg (or png) and renamed it. Everything is in ```QAA_report_plots```
    - copy the contents of the images directory back into personal computer, just for personal reference (and if I want to fix up the plot images on ```QAA_report.rmd``` and ```QAA_report.pdf```):
    
    ```
    cd ~/bioinfo/Bi623/PS/PS2

    scp -r joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/QAA_report_plots .
    ```

    - Just in case I do want to fix the way the plots look on ```QAA_report```: 
        - On my personal computer, find a way to use svgs in the Google Slides, and then make those Google Slides into pngs again, just like I did earlier, and replace those quality plots pngs in the ```QAA_report_plots``` folder. Update ```QAA_report.rmd``` and ```QAA_report.pdf```.
        - Then, upload ```QAA_report.pdf``` and the images folder ```QAA_report_plots``` from my personal computer onto Talapas, then onto Git:
        ```
        cd ~/bioinfo/Bi623/PS/PS2

        scp QAA_report* joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/

        # to copy the images directory and all its content:
        scp -r QAA_report_plots joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/
        ```

- counts files generated from htseq-count, copied from ```htseq_count``` directory into ```/projects/bgmp/joycew/bioinfo/Bi623/QAA/htseq_count_files```, and renamed properly: CHECK
- PDF report (turned into both GitHub, at a top level, and Canvas) CHECK

- Talapas batch script/code: CHECK
    - **git commit -m "copies of every single script used-- note that I made copies of scripts found in many different directories, so the paths to files might be wonky"**

- lab notebook: CHECK
    - Since this lab notebook has soooo much code, copy this ```Bi623_PS2_lab_notebook.md``` into Talapas, and then upload it.
    - Remember, this will be uploaded into the Bi621 lab notebook git repo too.

    ```
    # when completely done with lab notebook:
    cd /Users/joycewang/bioinfo/Bi621/Bi621_lab_notebook_joyce

    scp Bi623_PS2_lab_notebook* joycew@login2.talapas.uoregon.edu:/projects/bgmp/joycew/bioinfo/Bi623/QAA/

    # git add, commit, and push to the Bi621 lab notebook repo
    # git add, commit, and push to the QAA repo
    ```

