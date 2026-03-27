from print_substmat import print_substmat
from pathlib import Path


def test_print_substmat_header_half(capsys):
    """Test that the header row is formatted correctly for half format."""
    mat_file = Path(__file__).parent / "blosum62_half.mat"
    print_substmat(mat_file, "half")
    captured = capsys.readouterr()
    # Header with 5 spaces at start, then 4 spaces between amino acids
    assert captured.out.startswith(
        "     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n"
    )


def test_print_substmat_data_row_half(capsys):
    """Test that a data row is formatted correctly with right-justified numbers."""
    mat_file = Path(__file__).parent / "blosum62_half.mat"
    print_substmat(mat_file, "half")
    captured = capsys.readouterr()
    # Check that numbers are right-justified with 3 spaces separator
    assert "E   -1    0    0    2   -4    2    5" in captured.out


def test_print_substmat_header_full(capsys):
    """Test that the header row is formatted correctly for full format."""
    mat_file = Path(__file__).parent / "blosum62_full.mat"
    print_substmat(mat_file, "full")
    captured = capsys.readouterr()
    # Header with 5 spaces at start, then 4 spaces between amino acids
    assert captured.out.startswith(
        "     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n"
    )


def test_print_substmat_data_row_full(capsys):
    """Test that a data row is formatted correctly with right-justified numbers for full format."""
    mat_file = Path(__file__).parent / "blosum62_full.mat"
    print_substmat(mat_file, "full")
    captured = capsys.readouterr()
    # Check that numbers are right-justified with 3 spaces separator
    assert "E   -1    0    0    2   -4    2    5" in captured.out


def test_print_substmat_default_format(capsys):
    """Test that the default format is 'full'."""
    mat_file = Path(__file__).parent / "blosum62_full.mat"
    print_substmat(mat_file)  # Should use 'full' by default
    captured = capsys.readouterr()
    assert captured.out.startswith(
        "     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n"
    )
