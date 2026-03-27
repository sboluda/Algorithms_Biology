"""
Considerations:
- sys.stdin to take line by line the input fastq file
- readline().strip()
- take line by line until we reach the forth line with the quality
- consider the 4 lines a fastq has
- to get the correct score --> ord(ASCII_code) - 33
"""
import sys

nreads = 0
rname = ""
while True:
    rname = sys.stdin.readline().strip()
    if rname == "": break
    seq = sys.stdin.readline().strip()
    plus = sys.stdin.readline().strip()
    if plus!= "+":
        print("malformed FASTQ", file = sys.stderr)
        sys.exit(1)
    # quality string
    qs = sys.stdin.readline().strip()
    qv = [ord(e) - 33 for e in qs]
    # average quality
    aq = sum(qv)/len(qv)
    print("{}\t{:.3f}".format(rname, aq))
    nreads += 1
print("{} read(s)".format(nreads))
# Printed in the stderror not the stdint