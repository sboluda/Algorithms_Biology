"""
Program: decode-flag.py

Description:
This program reads an integer value from the command line representing
a SAM FLAG and prints the corresponding FLAG names.

Usage:
    python decode-flag.py <flag_value>

Example:
    python decode-flag.py 272

Output:
    REVERSE
    SECONDARY

The program should interpret the input integer as a bitwise FLAG and
print the names of all flags that are set in the given value, each on
a separate line.
"""

import sys

names = ['PAIRED', 'PROPER_PAIR', 'UNMAPPED', 'MATE_UNMAPPED', 'REVERSE', 'MATE_REVERSE', 'FIRST_IN_PAIR', 'SECOND_IN_PAIR', 'SECONDARY', 'QC_FAIL', 'DUPLICATE', 'SUPPLEMENTARY']

flag = int(sys.argv[1])
for i in range(12):
    if flag & (2**i) == 2**i: # We play with bits directly: If _ AND 1 == 1 --> It means that flag is active, because it must be 1 & 1, NOT 0 & 1
        print(names[i])


# flag = sys.argv[1]
# def decode_flag(flag):
#     """
#     Given a SAM FLAG integer, return a list called names with active flags in uppercase.

#     >>> decode_flag(99)
#     ['PAIRED', 'PROPER_PAIR', 'MATE_REVERSE', 'FIRST_IN_PAIR']

#     >>> decode_flag(147)
#     ['PAIRED', 'PROPER_PAIR', 'REVERSE', 'SECOND_IN_PAIR']

#     >>> decode_flag(4)
#     ['UNMAPPED']
#     """
#     names = ['PAIRED', 'PROPER_PAIR', 'UNMAPPED', 'MATE_UNMAPPED', 'REVERSE', 'MATE_REVERSE', 'FIRST_IN_PAIR', 'SECOND_IN_PAIR', 'SECONDARY', 'QC_FAIL', 'DUPLICATE', 'SUPPLEMENTARY']
#     flags = []
#     binary_values = bin(flag)
#     for i in range(12):
#         if binary_values[-1 - i] == "b":
#             break
#         elif binary_values[-1 - i] == "1":
#             flags.append(names[i])
#         else:
#             pass
#     return flags
