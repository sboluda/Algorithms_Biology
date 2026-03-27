def detect_format(file_path):
    """Detects alignment format based on file extension.

    >>> detect_format("test.fasta")
    'fasta'
    >>> detect_format("test.fa")
    'fasta'
    >>> detect_format("test.sto")
    'sto'
    >>> detect_format("test.aln")
    'aln'
    >>> detect_format("test.clustal")
    'aln'
    >>> detect_format("TEST.FASTA")
    'fasta'
    """
