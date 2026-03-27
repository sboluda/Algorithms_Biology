import sys
import os

# Constants
FASTA_LINE_WIDTH = 60  # Standard line width for FASTA and CLUSTALformat
ALN_ID_WIDTH = 30  # Width for sequence IDs in CLUSTAL and STOCKHOLM formats


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
    extension = file_path.strip().split(".")[-1].lower()
    if extension == "fa":
        return "fasta"
    elif extension == "clustal":
        return "aln"
    else:
        return extension


def read_alignment(file_path):
    """Reads an alignment file and returns a dictionary {seq_id: sequence}."""
    fmt = detect_format(file_path)
    if fmt == "fasta":
        return _read_fasta(file_path)
    if fmt == "aln":
        return _read_clustal(file_path)
    else:
        return _read_stockholm(file_path)

def _read_fasta(f):
    """Read FASTA format efficiently using lists."""
    sequences = {}
    seq_id = None
    seq_parts = []
    with open(f, "rt") as fasta:
        # We capture the first header
        seq_id = fasta.readline().strip()
        if not seq_id.startswith(">"):
            raise ValueError("Invalid FASTA format: first line must start with '>'")
        
        for line in fasta:
            if line.startswith(">"):                        # Start of a sequence
                sequences[seq_id] = "".join(seq_parts)
                seq_parts = []
                seq_id = line.strip()
            else:                                           # Sequence line
                seq_parts.append(line.strip())
    if seq_parts != []:                                     # Store the last sequence
        sequences[seq_id] = "".join(seq_parts)

    return sequences

def _read_stockholm(f):
    """Read Stockholm format efficiently using lists."""
    sequences = {}
    seq_id = None
    with open(f, "rt") as aln:
        for line in aln:
            if not line.startswith("#") and not line.startswith("//"): # Sequence line
                parts = line.strip().split()
                seq_id = parts[0]
                if seq_id not in sequences:
                    sequences[seq_id] = [parts[1]]
                else:
                    sequences[seq_id].append(parts[1])
    for keys, values in sequences.items():
        sequences[keys] = "".join(values)

    return sequences


def _read_clustal(f):
    """Read Clustal format efficiently using lists."""
    sequences = {}
    seq_id = None
    with open(f, "rt") as clustal:
        clustal.readline() # Skip the first line (header)
        for line in clustal:
            if not line.startswith(" "): # Ignore onsensus lines
                parts = line.strip().split()
                seq_id = parts[0]
                if seq_id not in sequences:
                    sequences[seq_id] = [parts[1]]
                else:
                    sequences[seq_id].append(parts[1])
    for keys, values in sequences.items():
        sequences[keys] = "".join(values)

    return sequences


def write_alignment(sequences, file_path):
    """Writes sequences to file in the detected format."""
    fmt = detect_format(file_path)
    if fmt == "fasta":
        return _write_fasta(file_path)
    if fmt == "aln":
        return _write_clustal(file_path)
    else:
        return _write_stockholm(file_path)

def _write_fasta(sequences, out):
    """Write FASTA format."""
    with open(out, "wt") as fasta:
        for key, sequence in sequences.items():
            # Header line (assuming it does not have a >)
            fasta.write(f'>{key}\n')

            # Sequence line (consider the length correction)
            for i in range(0, len(sequence), FASTA_LINE_WIDTH):
                fasta.write(sequence[i:i+FASTA_LINE_WIDTH] + "\n")

def _write_stockholm(sequences, out):
    """Write Stockholm format."""
    with open(out, "wt") as aln:
        # Write header of file line
        aln.write("# STOCKHOLM 1.0\n\n")
        seq_length = max(len(seq) for seq in sequences.values())
        for i in range(0, seq_length, ALN_ID_WIDTH):
            for key, sequence in sequences.items():
                aln.write(f'{key}\t{sequence[i:i+ALN_ID_WIDTH]}\n')
            aln.write("\n") # New line after each block
        aln.write("//\n")

def _write_clustal(sequences, out):
    """Write Clustal format."""
    # missing code
