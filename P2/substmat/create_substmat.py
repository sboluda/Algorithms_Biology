from pathlib import Path
from collections import Counter
import math

def create_substmat(filename, alphabet=None):
    counting_pairs = Counter()
    with open(filename, "rt") as fasta:
        if alphabet is None: # In case they do not provide an alphabet
            print("The alignment contains amino acids that are not in the provided alphabet")
            return None
        
        sequences = []
        for line in fasta:
            if line.startswith(">"):
                pass
            else:
                sequences.append(line.strip())
        
        mat_columns = [["_" for _ in range(len(sequences))] for _ in  range(len(sequences[0]))]
        
        count_bases = Counter() # To later count the total number of bases in all sequences
        for i in range(len(sequences)): # For loop to get ready the sequences in rows for later combinations computation
            for j in range(len(sequences[i])):
                mat_columns[j][i] = sequences[i][j]
                count_bases[sequences[i][j]] += 1

        for base in count_bases: # In case the alignment contains aa not present in the alphabet
            if base not in alphabet:
                print("The alignment contains amino acids that are not in the provided alphabet")
                return None
        

        count_pairs = Counter() # To later count the total number of pairs
        for combination in mat_columns: # Computation of all the combinations storing them in a dictionary
            for i in range(len(combination)-1):
                for j in range(i+1, len(combination)):
                    pair = combination[i] + combination[j]
                    count_pairs["".join(sorted(pair))] += 1
                
        total_pairs = sum(x for x in count_pairs.values()) # To compute the total number of all pairs
        observed_frequencies = {base: (value/total_pairs) for base, value in count_pairs.items()}


        total_bases = sum(x for x in count_bases.values()) # To compute the total number of all bases
        expected_frequences_single = {x: (y/total_bases) for x, y in count_bases.items()}
        
        expected_frequences_pairs = {}
        for pair in count_pairs:
            frequency = 0
            if pair != pair[::-1]:
                frequency = 2 * expected_frequences_single[pair[0]] * expected_frequences_single[pair[1]]
            else:
                frequency = expected_frequences_single[pair[0]]**2

            expected_frequences_pairs[pair] = frequency

        odds_ratio_10 = {pair: int(round(math.log10(observed_frequencies[pair] / expected_frequences_pairs[pair]) * 10, 0)) for pair in count_pairs}
        
        # Header
        print("     " + "    ".join(alphabet))

        # For each row in alphabet
        for i, res1 in enumerate(alphabet):
            row_values = res1
            for j, res2 in enumerate(alphabet[: i + 1]):  # only lower-triangular + diagonal
                pair = "".join(sorted([res1, res2]))
                mat_value = odds_ratio_10[pair]
                if mat_value < 0:
                    row_values += "   " + str(mat_value)
                else:
                    row_values += "    " + str(mat_value)
            print(row_values)
        
    return None
