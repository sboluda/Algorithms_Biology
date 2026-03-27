"""
Program: gap-compressed-identity.py

Description:
This program computes the gap-compressed identity of an alignment based on a CIGAR string.

Definition of gap-compressed identity:
- Consecutive indels (insertions or deletions) are counted as a single gap event
- Each gap event (regardless of its length) counts as ONE difference
- Matches (M) contribute to the number of matches
- Mismatches (X) count as mismatches
- Insertions (I) and deletions (D) contribute to gap events

Identity is computed as:
    Identity = matches / (matches + mismatches + gap_events)

Usage:
    python gap-compressed-identity.py <CIGAR_string>

Example:
    python gap-compressed-identity.py 3M2X4M2D2M3I

Output:
    matches: 9
    mismatches: 2
    openings: 2
    id: 0.69

The program parses the CIGAR string and groups consecutive indels (I or D) into a single gap event.
"""
import sys

cigar = sys.argv[1]
if len(cigar) % 2 != 0: # CIGAR must be even length
    print("Error: CIGAR must be an even length string", file=sys.stderr)
    sys.exit(1)

matches, mismatches, gaps = 0, 0, 0
for i in range(0, len(cigar), 2):
    # We assume the CIGAR can contain other type of operations (which we do not care about)
    event = cigar[i + 1]
    value = int(cigar[i])

    if event == "M": # Match
        matches += value
    elif event == "X": # Mismatch
        mismatches += value
    elif event == "I" or event == "D": # Gap (Insertion or Deletion)
        gaps += 1
    
id = round((matches / (matches + mismatches + gaps)), 2)

print(f"matches: {matches}")
print(f"mismatches: {mismatches}")
print(f"openings: {gaps}")
print(f"id: {id}", end="")