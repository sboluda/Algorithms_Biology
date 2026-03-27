def kmers(sequence, length):
    """
    Write a function that returns a list with all the words of length k contained in
    the query.

    >>> kmers("QLNFQLMSAGQLQ", 3)
    ['QLN', 'LNF', 'NFQ', 'FQL', 'QLM', 'LMS', 'MSA', 'SAG', 'AGQ', 'GQL', 'QLQ']
    >>> kmers("QLNFQLMSAGQLQ", 4)
    ['QLNF', 'LNFQ', 'NFQL', 'FQLM', 'QLMS', 'LMSA', 'MSAG', 'SAGQ', 'AGQL', 'GQLQ']"""
    ls_kmers = []
    for i in range(len(sequence) - length + 1):
        ls_kmers.append(sequence[i: i+length])
    return ls_kmers