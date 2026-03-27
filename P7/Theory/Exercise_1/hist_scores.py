# Use: $ zcat TESTX_H7YRLADXX_S1_L001_R1_001.fastq.gz | python avg-qual-self.py | python hist_scores.py
# Double pipe | to pass the data from script to script

import matplotlib.pyplot as plt
import sys

scores = []
# with open ("avg_scores.txt","rt") as fp:
#     for line in fp:
#         values = line.strip().split()
#         if len(values) == 2:
#             break
#         scores.append(float(values[2]))
for line in sys.stdin:
    read = line.strip().split("\t")
    scores.append(float(read[1]))

fig, ax = plt.subplots(1, 1)                # We create the fig and ax
ax.hist(scores, color = "#3894d2",        # We create the hist with the list scores
        edgecolor="white",
        linewidth=1.2)
plt.show()