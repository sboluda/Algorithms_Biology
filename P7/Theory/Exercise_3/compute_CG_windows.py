"""
BED format: 
    - Describes regions in the genome - 3 columns: Chr (contig), start, end
    - 0-based = Half-open coordinates --> [2, 4) Start included, end NOT
    - 1-bases = Closed coordinates --> [2, 4] 

Commands used:
gunzip <genome.fna.gz> > genome.fna.bgz
samtools faidx genome.fna.bgz
"""

import sys 
import matplotlib.pyplot as plt 

cg_reads = [] 

with open(sys.argv[1], "rt") as fasta: 
    # Skip the first header line 
    first_line = fasta.readline().strip() 
    if not first_line.startswith(">"): 
        print("Error: Fasta file does NOT start with a header", file=sys.stderr) 
        sys.exit(1)
    
    gc_count = 0 
    read_length = 0 
    for line in fasta: 
        line = line.strip().upper() 
        if line.startswith(">"): 
            if read_length > 0: # Avoid possible read_length == 0 situations 
                cg_reads.append(round((gc_count/read_length), 3)) 
            gc_count = 0 
            read_length = 0 
        
        else: 
            read_length += len(line) 
            gc_count += sum( 
                1 for nucleotide in 
                line if nucleotide in "CG" 
                ) 
    if read_length > 0:
        cg_reads.append(round((gc_count/read_length), 2))  # Compute the last read 
        
fig, ax = plt.subplots(1, 1) # We create the fig and ax 
ax.hist(cg_reads, color = "#3894d2", # We create the hist with the list scores 
        edgecolor="white", 
        linewidth=1.2) 
plt.show()