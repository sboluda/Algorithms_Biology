"""
Program: parse-meth.py

Description:
This program parses methylation information (MM/ML tags) and prints the positions of methylated cytosines along with their corresponding methylation likelihoods.

Input:
- A DNA sequence
- An MM string (comma-separated counts of skipped positions)
- An ML string (comma-separated likelihood values)

All inputs are provided as command-line arguments.

Usage:
    python parse-meth.py <sequence> <MM> <ML>

Example:
    python parse-meth.py ACGTTAACGTTCCCGATCT 0,0,2 200,215,199

Output:
    Prints one line per methylated cytosine in the format:
        position:likelihood

    Example:
        1:200
        7:215
        13:199

The program maps the MM counts to positions in the sequence and associates each methylation event with its corresponding likelihood from the ML string.
"""
import sys


sequence = sys.argv[1].strip()                         # We store the sequence
mm = list(map(int, sys.argv[2].strip().split(",")))    # We store the MM (C methylated positions) - Use map() to apply int() to all elements in the list and list() to create the list
ml = list(map(int, sys.argv[3].strip().split(",")))    # We store the ML (methylation likelihood of each C)
mm_length = len(mm)
m_index = 0

for i in range(len(sequence)):
    if sequence[i] == "C":
        if m_index == mm_length:          # No more methylated Cs
            break

        elif mm[m_index] == 0:              # It's a match
            print(f"{i}:{ml[m_index]}")
            m_index += 1
        else:                             # Substract 1 until reaching 0 (the methylated C)
            mm[m_index] -= 1



