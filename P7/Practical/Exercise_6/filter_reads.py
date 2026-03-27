"""
Program: filter_reads.py

Description:
This program filters a SAM file, keeping only reads that satisfy the following conditions:
- Mapping quality is above a given threshold (> 30)
- The read is NOT marked as SECONDARY
- SAM header lines are preserved and printed without filtering

Usage:
    cat input.sam | python filter_reads.py <quality_threshold> | wc -l

Input:
    - SAM file read from standard input (stdin)
    - Quality threshold provided as a command-line argument

Output:
    - Filtered SAM records printed to standard output
    - Header lines (starting with "@") are always printed
    - Only reads with:
        * mapping quality > threshold
        * not SECONDARY (FLAG 256 not set)
      are included in the output

### With the previously computed SAM file, the output should be: 287 ###
"""

import sys

thr = int(sys.argv[1])

for line in sys.stdin:
    line = line.strip()
    if line.startswith("@"):
        print(line)
    else:
        read = line.split()
        flag = int(read[1])
        q_score = int(read[4])

        if (flag & 256) == 0 and q_score > thr:
            print(line)

# Teacher solution:

'''
thr = int(sys.argv[1])

for line in sys.stdin:
    line = line.strip()

    if line[0] == "@":
        print(line)
        continue
        
    fields = line.split("\t")
    rname = fields[0]
    flags = int(fields[1])
    qual = int(fields[4])

    if (flags & 256) != 0:
        # print("{} SECONDARY".format(rname), file=sys.stderr)
        continue
    if qual <= thr:
        continue
    print(line)
'''
