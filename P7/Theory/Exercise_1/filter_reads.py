"""
1. Filter paired-end FASTQ files based on the average quality score of each read.
2. Two input FASTQ files and two output FASTQ files as command-line arguments.
3. Write reads to the output files only if both reads in a pair have an average quality score of 20 or greater.

Usage: python filter-fastq.py <input1.fastq> <input2.fastq> <output1.fastq> <output2.fastq>

To later gz the files -> gzip <file_name>
"""
import argparse

parser = argparse.ArgumentParser(
    description = "Given two paired end fastq files - Write reads to new fastq files only if both reads in a pair have an average quality score of 20 or greater."
)

parser.add_argument("first_fastq_input", type=str, help="First pair of the fastq (R1)")
parser.add_argument("second_fastq_input", type=str, help="Second pair of the fastq (R2)")
parser.add_argument("first_fastq_output", type=str, help="First pair of the fastq after filter")
parser.add_argument("second_fastq_output", type=str, help="Second pair of the fastq after filter")

args = parser.parse_args()

import sys
import gzip # Alternative to read zip files without using zcat in the shell

first_input_fastq = args.first_fastq_input
second_input_fastq = args.second_fastq_input
first_output_fastq = args.first_fastq_output
second_output_fastq = args.second_fastq_output
total_reads = 0

with (gzip.open(first_input_fastq, "rt") as input1,
      gzip.open(second_input_fastq, "rt") as input2,
           open(first_output_fastq, "wt") as output1,
           open(second_output_fastq, "wt") as output2): # We can use with open of multiple files!
    
    count = 1
    for line1, line2 in zip(input1, input2):
        if count%4 == 1:     # Header line
            header_line_1 = line1   # We do NOT use strip() to conserve the \n
            header_line_2 = line2   # We do NOT use strip() to conserve the \n
        
        elif count%4 == 2:   # Sequence line
            seq_line_1 = line1  # We do NOT use strip() to conserve the \n
            seq_line_2 = line2  # We do NOT use strip() to conserve the \n
            if len(seq_line_1.strip()) != len(seq_line_2.strip()):
                print("ERROR: Pair-end reads are not of the same length!")
                sys.exit(1)

        elif count%4 == 3:   # + line
            plus_line_1 = line1
            plus_line_2 = line2

        elif count%4 == 0:   # Quality line
            quality_line_1 = line1
            quality_line_2 = line2

            if len(seq_line_1.strip()) != len(quality_line_1.strip()):
                print("ERROR: Sequence and quality length lines are different")
            
            average_phred_1 = sum(ord(chr) - 33 for chr in quality_line_1.strip()) / len(quality_line_1.strip())
            average_phred_2 = sum(ord(chr) - 33 for chr in quality_line_2.strip()) / len(quality_line_2.strip())
            
            if average_phred_1 < 20 or average_phred_2 < 20:
                pass
            else: # They pass the threshold, so we put them in the new fastq
                output1.write(header_line_1 +
                              seq_line_1 + 
                              plus_line_1 +
                              quality_line_1)

                output2.write(header_line_2 +
                              seq_line_2 + 
                              plus_line_2 +
                              quality_line_2)
                total_reads += 1   
        count += 1

print(f"TOTAL READS: {total_reads}")