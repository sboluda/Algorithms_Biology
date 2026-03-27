"""
Calculate and visualize the average Phred quality score per
position across multiple FASTQ reads.
1. Read FASTQ data from standard input
2. Aggregate quality scores for each base position
3. Compute the average quality score for each position, print these averages, and generate a plot of these averages.
4. Usage: zcat reads.fastq.gz python avg-qual-by-pos.py
"""
import sys
import matplotlib.pyplot as plt

count = 1
mapping_quality = {}

for line in sys.stdin:
    if count%4 == 0: # Quality line
        quality_sequence = line.strip()
        for position in range(len(quality_sequence)): # We iterate over all individual bases
            if position not in mapping_quality:
                mapping_quality[position] = []
            mapping_quality[position].append(ord(quality_sequence[position]) - 33)
    count += 1

# Once finished, we have all scores for each position, now we compute the average
avg_phred_position = []
positions = []

for key in sorted(mapping_quality.keys()):
    positions.append(key)
    list_values = mapping_quality[key]
    avg_phred_position.append(round(sum(list_values) / len(list_values), 3))

fig, ax = plt.subplots(1, 1)                # We create the fig and ax
ax.plot(positions, avg_phred_position, color = "#3894d2")
ax.set_ylim(25.0, 45.0)

plt.show()
