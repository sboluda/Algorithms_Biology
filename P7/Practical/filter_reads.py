# Print reads with quality > 30 and NOT secondary
# Usage: 
# cat file.sam | python filter_reads.py qualty_value | wc -l

import sys

for line in sys.stdin:
    line = line.strip()
    if line.startswith("@"):
        print(line)
    else:
        read = line.split()
        flag = int(read[1])
        q_score = int(read[4])
        if (flag & 256 == 0) and q_score > int(sys.argv[1]):
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
    if (flags & 256) == 256:
        print("{} SECONDARY".format(rname), file=sys.stderr)
        continue
    if qual > thr:
        continue
    print(line)
'''
