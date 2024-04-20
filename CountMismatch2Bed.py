#!/usr/bin/env python
import subprocess as sp
import shlex
import pysam
import sys
import re

####Summary: This script was used to call mismatches from bam file for DMS-MaPseq. Mismatches happened 3 nt within indels were discarded.
####         The output is a bed file with reference nucleotide and mismatch number for each mismatched location
####Author: Changhao Li
####E-mail: chli.bioinfo@gmail.com
####Reference: based on mismatch_counting.py written by Meiyue Wang from 2019 Methods paper

bamfile = sys.argv[1]

baseInfo = {}
for read in pysam.AlignmentFile(bamfile, 'r'):
    readname    = read.query_name
    chrom       = read.reference_name
    pos         = int(read.pos) + 1
    MD          = read.get_tag("MD")
    md_match    = re.findall(r'\d+',MD)
    md_mismatch = re.findall(r'\D+',MD)
    cigar       = read.cigarstring
    cg_num      = re.findall(r'\d+',cigar)
    cg_lt       = re.findall(r'\D+',cigar)



    region_range_type = {}  #### used to store [start, end, region_type] and the key is number
    indel_reg         = []  #### used to store up/down 3 nt of indel regions
    crt_pos = pos #crt represents current
    #### matched region, indel up/downstream 3 bp region and intron region
    for n in range(len(cg_lt)):
        if cg_lt[n] == "M":    #### for matched regions
            region_range_type[n] = [crt_pos,crt_pos+int(cg_num[n])-1,"M"]
            crt_pos += int(cg_num[n])
        elif cg_lt[n] == "D":  #### for indel regions; For deletion, there may be multiple deletions, so we should store the start and end of the deletion
            region_range_type[n] = [crt_pos, crt_pos+int(cg_num[n])-1, "ID"] #### ID means indel
            indel_reg += range(crt_pos - 3, crt_pos + int(cg_num[n]) + 4)
            crt_pos += int(cg_num[n])
        elif cg_lt[n] == "I":  #### for indel regions; Different from deletion, insertion only happened in one locate, so we store just one postion 
            region_range_type[n] = [crt_pos, crt_pos,"ID"]
            indel_reg += range(crt_pos - 3, crt_pos + 4)
        elif cg_lt[n] == "N":  #### for intron regions or skipped regions
            region_range_type[n] = [crt_pos, crt_pos+int(cg_num[n])-1, "N"]
            crt_pos += int(cg_num[n])
    

    ##### get mismatch for each position
    crt_pos = pos
    region_initiation = 0
    for n in range(len(md_mismatch)):

        ref = md_mismatch[n]  #### mismatch happened reference nucleotide
        base = crt_pos + int(md_match[n]) #### base: mismatch happened position
        region_i = region_initiation      #### region_initiation can be changed when move to the next region (The else part and the if part with region type N)

        while region_i < len(region_range_type):    #### iterate each cigar regions
            if(base >= region_range_type[region_i][0] and base <= region_range_type[region_i][1]):  #### The if part, which mean the mismatch happened in this region
                if(region_range_type[region_i][2] == "M"):
                    crt_pos = base + len(re.sub("\^", "", md_mismatch[n]))
                    ####count this base
                    tmpkey = "*".join([chrom,str(base)])
                    if base not in indel_reg:
                        try:
                            baseInfo[tmpkey][1] +=1
                        except KeyError:
                            baseInfo[tmpkey] = [ref, 1]

                    region_initiation = region_i
                    break

                if(region_range_type[region_i][2] == "N"):
                    #### the base should not be in intron, we should relocate it 
                    rel_len = base-region_range_type[region_i][0]
                    base    = region_range_type[region_i][1] + rel_len + 1
                    crt_pos = base + len(re.sub("\^", "", md_mismatch[n]))

                    region_i+=1
                    region_initiation = region_i
                    continue
                
                if(region_range_type[region_i][2] == "ID"):
                    #### the base was in indel region

                    crt_pos = base + len(re.sub("\^", "", md_mismatch[n]))
                    region_initiation = region_i
                    break

            else: #### The else part, which means the mismatch happened downstream current region
                if(region_range_type[region_i][2] == "M"):

                    region_i+=1
                    region_initiation = region_i
                    continue
                if(region_range_type[region_i][2] == "N"):
                    rel_len = base-region_range_type[region_i][0]
                    base    = region_range_type[region_i][1] + rel_len + 1
                    crt_pos = base + len(re.sub("\^", "", md_mismatch[n]))

                    region_i+=1
                    region_initiation = region_i
                    continue
                
                if(region_range_type[region_i][2] == "ID"):

                    region_i+=1
                    region_initiation = region_i
                    continue
                
                continue
                    

        

for key in baseInfo:
    tmpkey = key.split("*")
    chrom = tmpkey[0]
    start = str(int(tmpkey[1]) - 1)
    end   = str(int(tmpkey[1]))
    name  = "_".join([chrom, end, baseInfo[key][0]])
    count = str(baseInfo[key][1])

    print chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + count
