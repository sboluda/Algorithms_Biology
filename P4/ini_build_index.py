def build_index(db_seqs, k):
    """
    Preprocess the database into a k-mer index.
    Keys are sorted alphabetically.

    Returns a dict:
        { kmer: [(indx, pos), ...], ... }

    >>> build_index(["MRTAY", "KTMKT"], 3)
    {'KTM': [(1, 0)], 'MKT': [(1, 2)], 'MRT': [(0, 0)], 'RTA': [(0, 1)], 'TAY': [(0, 2)], 'TMK': [(1, 1)]}
    >>> build_index(["AAA"], 2)
    {'AA': [(0, 0), (0, 1)]}
    """
    db = {}
    for i, seq in enumerate(db_seqs): # Using enumarete() allows us to obtain both the position of the seq and the seq itself
        for j in range(len(seq) - k + 1): # To have a "sliding window" frame
            # Compute the substring (remember to apply the window [j: j + k])
            substring = seq[j: j + k]
            
            # We compute the tuple, which is just (i, j)
            index = (i, j)

            if substring not in db:
                # Store the new kmer in the dictionary
                db[substring] = [index]
            else:
                db[substring].append(index)
    return dict(sorted(db.items()))

if __name__ == "__main__":
    print(build_index(["MRTAY", "KTMKT"], 3))
