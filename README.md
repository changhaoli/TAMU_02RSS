# CountMismatch2Bed.py 
CountMismatch2Bed.py is a python2 script, which can be used to call mismatches for DMS-MaPseq data.
CountMismatch2Bed.py only needs a sorted input.bam file and outputs a bed file that includes location and mismatch count information.

## Requirement
python2 >= 2.7.10\
pysam >= 0.9.1.4

## Usage
python CountMismatch2Bed.py input.bam > output.bed
