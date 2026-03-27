from Bio.Align import substitution_matrices

def score_words(words):
    """
    Write a function that takes a list of words (all words with the same length) and performs
    all-to-all comparisons of words using a BLOSUM62 substitution matrix. The output should
    be as indicated in the tests: a dictionary of with the words as keys and a list of tuples as
    values that consist of the score and the aligned word, sorted from larger to smaller score.

    >>> score_words(['AAA', 'AAS', 'ASA', 'SSA', 'SSS'])
    {'AAA': [(12.0, 'AAA'), (9.0, 'ASA'), (9.0, 'AAS'), (6.0, 'SSA'), (3.0, 'SSS')], 'AAS': [(12.0, 'AAS'), (9.0, 'AAA'), (6.0, 'SSS'), (6.0, 'ASA'), (3.0, 'SSA')], 'ASA': [(12.0, 'ASA'), (9.0, 'SSA'), (9.0, 'AAA'), (6.0, 'SSS'), (6.0, 'AAS')], 'SSA': [(12.0, 'SSA'), (9.0, 'SSS'), (9.0, 'ASA'), (6.0, 'AAA'), (3.0, 'AAS')], 'SSS': [(12.0, 'SSS'), (9.0, 'SSA'), (6.0, 'ASA'), (6.0, 'AAS'), (3.0, 'AAA')]}
    """
    subst_mat = substitution_matrices.load("BLOSUM62")

    scores = {}
    for main_kmer in words:
        comparisons = []
        for comparison_kmer in words:
            score = sum(
                subst_mat[i, j]
                for i, j in zip(main_kmer, comparison_kmer)
            )
            comparisons.append((score, comparison_kmer))
        scores[main_kmer] = sorted(comparisons, key=lambda x: (x[0], x[1]), reverse=True)
    return scores