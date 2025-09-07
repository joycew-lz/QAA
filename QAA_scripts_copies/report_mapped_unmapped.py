#!/usr/bin/env python

mapped_count = 0
unmapped_count = 0

with open("/STAR_Align_and_picard/SRR25630394_Aligned.sorted.dedup.sam", 'r') as f:
    for line in f:
        if line.startswith("@"):
            continue # skip next step when line is a header line
        
        alignment_line = line.strip('\n').split() # returns alignment line list for current read,
        # every item of the list is split by a space
        flag = int(alignment_line[1]) # get the flag in this line

        # for bit 256: check if secondary alignment flag is true to see if current read has multiple alignments;
        # do NOT count reads as aligned more than once.
        if((flag & 256) == 256):
            # bit 256 is true, so the read IS a secondary alignment
            continue # skip next steps when read IS a secondary alignment (when read has multiple alignments)
        # if flag & 256 != 256, then bit 256 is false, and we can count
        # if the current read is mapped to the reference genome

        # for bit 4: check if current read is mapped to the reference genome:
        if((flag & 4) != 4):
            # bit 4 is false, so mapped = True
            mapped_count += 1
        else:
            # bit 4 is true, so mapped = False
            unmapped_count += 1
    
    print("Number of MAPPED reads:", mapped_count)
    print("Number of UNMAPPED reads:", unmapped_count)
        
        
