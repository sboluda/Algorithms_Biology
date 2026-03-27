"""
Reservoir Sampling Algorithm (for FASTQ reads)

1. Initialize a list (reservoir) to store k selected reads.

2. For the first k reads encountered:
   - Add them directly to the reservoir.
   - These occupy indices 0 to k-1.

3. For each subsequent read (the i-th read, where i >= k):
   - Generate a random integer j between 0 and i-1 (inclusive).
   - If j < k:
       Replace the element at index j in the reservoir with the new read.
   - Otherwise:
       Discard the new read.

This algorithm ensures that each read from the input stream
has an equal probability of being included in the final sample,
while using only O(k) memory and a single pass over the data.

Usage:
zcat reads.fastq.gz | python select-random-reads.py <k>
"""

import sys
import random

k = int(sys.argv[1])

reads = []
count = 1
i = 0

for line in sys.stdin:
    if count % 4 == 1: # Header line
        if not line.strip().startswith("@"):
            print("Error: file does not have a correct fastq format", file=sys.stderr)
            sys.exit(1)
        read = line
    
    elif count % 4 == 2: # Sequence line
        read += line
    
    elif count % 4 == 3: # + line
        read += line

    else: # Quality line
        read += line

        if i >= k:
            j = random.randint(0, i-1)
            if j < k:
                reads[j] = read
        else:
            reads.append(read)
        i += 1
    count += 1

# We use += line in order to store the \n of the whole read

for i in range(len(reads)):
    print(f"Read {i + 1}:")
    print(reads[i])