from Bio import SeqIO


def count_matches_fasta(fasta_file):
    """
    >>> count_matches_fasta("fastcats.fasta")
    5
    """
    records = SeqIO.parse(fasta_file, "fasta")
    seqs = [record.seq for record in records]

    # FILL IN


if __name__ == "__main__":
    print(count_matches_fasta("fastcats.fasta"))
