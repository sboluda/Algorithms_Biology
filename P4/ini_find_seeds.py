from Bio.Align import substitution_matrices
from ini_build_index import build_index

# from build_index import build_index

subst_mat = substitution_matrices.load("BLOSUM62")


def find_seeds(query, db, k, t):
    """
    Return all seeds whose BLOSUM62 score against a query k-mer meets
    threshold t, using a prebuilt k-mer index for lookups.

    >>> db = {'KTM': [(1, 0)], 'MKT': [(1, 2)], 'MRT': [(0, 0)], 'RTA': [(0, 1)], 'TAY': [(0, 2)], 'TMK': [(1, 1)]}
    >>> find_seeds("MKT", db, 3, 11)
    [{'db_iseq': 1, 'db_pos': 2, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MKT', 'score': 15.0}, {'db_iseq': 0, 'db_pos': 0, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MRT', 'score': 12.0}]
    """
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

if __name__ == "__main__":
    k = 3
    t = 11
    db = {
        "KTM": [(1, 0)],
        "MKT": [(1, 2)],
        "MRT": [(0, 0)],
        "RTA": [(0, 1)],
        "TAY": [(0, 2)],
        "TMK": [(1, 1)],
    }
    print(find_seeds("MKT", db, k, t))
