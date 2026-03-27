from Bio import SeqIO


def count_matches_stars(fasta_file):
    """
    >>> count_matches_stars("fastcats.fasta")
    '       ** * **'
    """
    records = SeqIO.parse(fasta_file, "fasta")
    seqs = [record.seq for record in records]

    # FILL IN CODE HERE


if __name__ == "__main__":
    print(count_matches_stars("fastcats.fasta"))
