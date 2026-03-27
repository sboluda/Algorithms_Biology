import sys
from tree_nodes import Node, newick2nodes
from ini_fasta2dict import fasta2dict
from ini_align_profiles_names import align_profiles

def node2splits_align(i_node, nodes, id2seq, gap):  #  this line requires a change
    """
    Aligns the profiles of the two subtrees of the node i_node, and returns
    the list of names of the sequences in the subtree rooted at i_node.

    >>> node2splits_align(0, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr', 'hmgt_mouse']
    >>> node2splits_align(1, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr']
    >>> node2splits_align(2, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat']
    """
    if nodes[i_node].left == 0:
        return [nodes[i_node].name]
    
    left_list = node2splits_align(
        nodes[i_node].left,
        nodes,
        id2seq,
        gap
        )
    
    right_list = node2splits_align(
        nodes[i_node].right,
        nodes,
        id2seq,
        gap
        )
    align_profiles(id2seq, left_list, right_list, gap)
    
    return left_list + right_list

if __name__ == "__main__":
    nodes = newick2nodes(sys.argv[1])  # "hmgb.dnd"
    id2seq = fasta2dict(sys.argv[2])  # "hmgb.fasta"
    node2splits_align(0, nodes, id2seq, -4)
    for name, seq in id2seq.items():
        print(f">{name}\n{seq}")

"""
Codes imported for the msa function:
def newick2nodes(newick):
    # newick can be a file path or a newick string
    if os.path.isfile(newick):
        # read from file
        with open(newick) as f:
            tree = Phylo.read(f, "newick")
    else:
        # read from string
        tree = Phylo.read(StringIO(newick), "newick")
    return _tree2nodes(tree)

def fasta2dict(filename):
    records = SeqIO.parse(filename, "fasta")
    sequences = {
        record.id: record.seq
        for record in records
    }
    return sequences

def align_profiles(profile_i, profile_j, gap):
    n_i = len(next(iter(profile_i.values())))
    n_j = len(next(iter(profile_j.values())))

    # Initialize scores and traceback matrices
    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    # Initialize edges with gaps
    for i in range(1, n_i + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1  # up (gap in profile_j)

    for j in range(1, n_j + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1  # left (gap in profile_i)

    # Fill the matrices
    for i in range(1, n_i + 1): # A loop for each sequence in profile i 
        for j in range(1, n_j + 1): # A loop for each sequence in profile j
            # Average substitution score across all sequence pairs, including gap penalties
            score = 0
            n_pairs = 0  # count of valid pairs (non-gap)
            for seq_i in profile_i.values():
                for seq_j in profile_j.values():
                    aa_i = seq_i[i - 1] # One particular aa in one of the sequences from profile i
                    aa_j = seq_j[j - 1] # One particular aa in one of the sequences from profile j
                    if aa_i == "-" and aa_j == "-":  # We first check if both are gaps (we EXCLUDE this value)
                        pass
                    elif aa_i == "-" or aa_j == "-": # We check for gap and aa case
                        score += gap
                        n_pairs += 1 
                    else:                            # We check for the rest
                        score += subst_mat[aa_i.upper(), aa_j.upper()]
                        n_pairs += 1
                    # We add a pair for all cases EXCEPT the gap and gap!!
            
            # Now we can compute the average score
            score = score / n_pairs if n_pairs > 0 else 0  # avoid division by zero

            diag = scores[i - 1][j - 1] + score  # match/mismatch
            left = scores[i][j - 1] + gap  # insertion (gap in profile_i)
            up = scores[i - 1][j] + gap  # deletion (gap in profile_j)

            scores[i][j] = max(diag, left, up)

            if diag >= left and diag >= up:
                traceback[i][j] = 0  # diagonal
            elif left >= up:
                traceback[i][j] = -1  # left
            else:
                traceback[i][j] = 1  # up

    # Traceback
    i, j = n_i, n_j
    aln = {name: [] for name in {**profile_i, **profile_j}}

    while i > 0 or j > 0:
        if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            for name, seq in profile_i.items():
                aln[name].append(seq[i])
            for name, seq in profile_j.items():
                aln[name].append(seq[j])
        elif traceback[i][j] == -1:  # left (gap in profile_i)
            j -= 1
            for name in profile_i:
                aln[name].append("-")
            for name, seq in profile_j.items():
                aln[name].append(seq[j])
        else:  # up (gap in profile_j)
            i -= 1
            for name, seq in profile_i.items():
                aln[name].append(seq[i])
            for name in profile_j:
                aln[name].append("-")

    return round(scores[-1][-1], 2), {
        name: "".join(reversed(chars)) for name, chars in aln.items()
    }
"""