from read_substmat_full import read_substmat_full
from read_substmat_half import read_substmat_half


def print_substmat(filename, format_type):
    """
    Prints a substitution matrix.

    Args:
        filename: Path to the substitution matrix file
        format_type: "full" or "half" to specify the matrix format
    """
    # YOUR CODE HERE


if __name__ == "__main__":
    print_substmat("blosum62_full.mat", "full")
    print("\n")
    print_substmat("blosum62_half.mat", "half")
