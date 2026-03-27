"""
Considerations:
- sys.stdin to take line by line the input fastq file
- readline().strip()
- take line by line until we reach the forth line with the quality
- consider the 4 lines a fastq has
- to get the correct score --> ord(ASCII_code) - 33

EXERCISE
1. Parse the 2 FASTQ files:
    TESTX_H7YRLADXX_S1_L001_R1_001.fastq.gz and
    TESTX_H7YRLADXX_S1_L001_R2_001.fastq.gz

2. Check that qualities have the same length than the sequences.
3. Compute mean quality of each read.
4. Read FASTQ formatted data from stdin
5. Calculate the average Phred quality score for each read,
6. Print the read name and average quality score to stdout

Usage:
zcat TESTX_H7YRLADXX_S1_L001_R1_001.fastq.gz | avg-qual-self.py
"""
import sys

count = 1
num_sequences = 0

for line in sys.stdin: # We collect all lines from the previous pipe
    if count%4 == 1:     # Header line
        new_line = line.strip() + "\t"
        num_sequences += 1
    
    elif count%4 == 2:   # Sequence line
        seq_length = len(line.strip())

    elif count%4 == 3:   # + line
        if line.strip() != "+":
            print("ERROR: No '+' line")
            sys.exit(1)

    elif count%4 == 0:
        seq_values = line.strip()
        if seq_length != len(seq_values):
            print("ERROR: Sequence and quality length lines are different")
        
        average = sum(ord(chr) - 33 for chr in seq_values) / len(seq_values)
        new_line += f'{average:.3f}'
        print(new_line) # Print to stdout
    count += 1    
# print(f'TOTAL: {num_sequences} sequences')