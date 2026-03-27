from Bio import SeqIO
from Bio.Align import substitution_matrices
import sys

subst_mat = substitution_matrices.load("BLOSUM62")

from ini_build_index import build_index
from ini_find_seeds import find_seeds
from ini_extend_seeds import extend_seeds
from ini_merge_overlapping import merge_overlapping

def build_index(db_seqs, k):
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

def find_seeds(query, db, k, t):
    subst_mat = substitution_matrices.load("BLOSUM62")
    query_kmers = db
    seeds = []
    for i in range(len(query) - k + 1):
        for key_kmer, positions in query_kmers.items():
            query_kmer = query[i: i + k]
            score = sum(
                subst_mat[i, j]
                for i, j in zip(query_kmer, key_kmer)
            )
            if score > t:
                for db_iseq, db_pos in positions:
                    seeds.append(
                        {
                            "db_iseq": db_iseq,
                            "db_pos": db_pos,
                            "q_pos": i,
                            "q_kmer": query_kmer,  # original query k-mer
                            "db_kmer": key_kmer,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds

def extend_seeds(seeds, query, db_seqs, k, X=15):
    hits = []  # will contain actual hits
    seen = (
        set()
    )  # to avoid duplicates, stores tuples (db_iseq, q_start, q_end, db_start, db_end)

    for s in seeds:
        q_pos = s["q_pos"]            # Query kmer match start position
        db_iseq = s["db_iseq"]        # Sequence position in DB (our db_seqs list)
        db_pos = s["db_pos"]          # Sequence kmer match start position
        db_seq = db_seqs[db_iseq]   # The Sequence itself

        # --- Extend right ---
        score, best_score, best_right = 0, 0, 0
        right = 0
        while (
            q_pos + k + right < len(query)
            and db_pos + k + right < len(db_seq)
            and score >= best_score - X
        ):
            score += subst_mat[                         # We compute the next right index position value for query and db_seq
                        query[q_pos + k + right] ,
                        db_seq[db_pos + k + right]
                        ] 
            if score > best_score:                      # If the score passes the best_score... 
                best_score = score                      # We store the new best_score
                best_right = right + 1                  # We store the new best_right position score
            right += 1                                  # We move one further position to the right

        # --- Extend left ---
        score, best_score, best_left = 0, 0, 0
        left = 0
        while (
            q_pos - left - 1 >= 0 
            and db_pos - left - 1 >= 0 
            and score >= best_score - X
        ):
            score += subst_mat[ 
                        query[q_pos - left - 1] ,
                        db_seq[db_pos - left - 1]
                        ] 
            if score > best_score:
                best_score = score
                best_left = left + 1
            left += 1

        q_start = q_pos - best_left
        q_end = q_pos + k + best_right - 1
        db_start = db_pos - best_left
        db_end = db_pos + k + best_right - 1

        key = (db_iseq, q_start, q_end, db_start, db_end)
        if key not in seen:
            seen.add(key)
            hits.append(
                {
                    "db_iseq": db_iseq,
                    "query_start": q_start,
                    "query_end": q_end,
                    "db_start": db_start,
                    "db_end": db_end,
                    "query_seq": query[q_start : q_end + 1],
                    "db_seq": db_seq[db_start : db_end + 1],
                    "seed": s["q_kmer"],
                }
            )

    return hits

def merge_overlapping(hits, query, db_seqs):
    merged_hits = []

    groups = {}
    for h in hits:
        diagonal = h["db_start"] - h["query_start"]  # To know the diagonal we must substract both db and query starts to reach the beginning of the diag 
        key = (h["db_iseq"], diagonal)
        if key not in groups:
            groups[key] = []
        groups[key].append(h)
    # We end up having a dictionary of hits grouped by:
        # If they belong to the same sequence from the DB
        # If they match the query in the same diagonal
    # We store the hits with these characteristics in common in the same key
    # We then store all the possible query combinations for all the diagonal and DB_sequences in our DB

    for (db_iseq, diagonal), group in groups.items():  # iterate over each db sequence and diagonal
        group.sort(
            key=lambda x: x["query_start"]      
        )  # sort hits by query_start to ensure they are in the correct order for merging
           # If we had (0, 0): [h1, h3, h5] --> [h3, h1, h5]

        current = dict(
            group[0]
        )  # start with the first hit in the group as the current hit to merge with subsequent hits
           # In the example from before, we store h3
           # We use dict() to create a COPY of h3 (but remember that h3 is already a dictionary)
        
        for h in group[1:]:
            if h["query_start"] <= current["query_end"] + 1:
                current["query_end"] = max(current["query_end"], h["query_end"])
                current["db_end"] = max(current["db_end"], h["db_end"])
                current["query_seq"] = query[
                    current["query_start"] : current["query_end"] + 1
                ]
                current["db_seq"] = db_seqs[db_iseq][
                    current["db_start"] : current["db_end"] + 1
                ]
            else:
                merged_hits.append(current) # store the hit
                current = dict(h) # start a new current hit with the next hit in the group
                                  # So we make a COPY of the hit we just compared, that is, comparing h1 with h5 

        merged_hits.append(current)

    return merged_hits

if __name__ == "__main__":
    if len(sys.argv) != 3: # Basic error handling
        print("Usage: python blast.py query.fasta database.fasta")
        sys.exit(1)

    k = 3  # k-mer length
    t = 11  # BLOSUM62 score threshold

    # We use the sys.argv[] to access the name of the file
    query_file = sys.argv[1]      
    db_file = sys.argv[2]

    # We use SeqIO to extract the sequences:
        # SeqIO.read for single sequences
        # .seq returns a Seq object -> Using str() we transform it into a string
        # In SeqIO.parse it's important to only focus on the .seq (that's why it's str(r.seq) and not just r)
    query = str(SeqIO.read(query_file, "fasta").seq)

    db_seqs = [str(r.seq) for r in SeqIO.parse(db_file, "fasta")]

    print(f"\nPreprocessing database (k={k})...")
    index = build_index(db_seqs, k)

    print(f"\nFinding seeds (threshold t={t})...")
    seeds = find_seeds(query, index, k, t)

    print(f"\nExtending seeds...")
    hits = extend_seeds(seeds, query, db_seqs, k)

    print(f"\nMerging overlapping hits...")
    merged = merge_overlapping(hits, query, db_seqs)

    print(f"\nMerged hits sorted by length (longest first): {len(merged)}")
    merged.sort(key=lambda h: h["query_end"] - h["query_start"], reverse=True)
    for h in merged[:10]:  # print top 10 longest hits
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_iseq']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")
