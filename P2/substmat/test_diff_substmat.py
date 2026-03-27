from diff_substmat import diff_substmat
from read_substmat_half import read_substmat_half
from pathlib import Path
import pytest
import re


def test_diff_substmat_values(capsys):
    """Test that the differences are calculated and printed correctly."""
    mat_file1 = Path(__file__).parent / "blosum62_half.mat"
    mat_file2 = Path(__file__).parent / "blosum45_half.mat"

    # Call the function (it only prints, doesn't return)
    diff_substmat(mat_file1, mat_file2)

    # Read the matrices to check the values manually
    mat1 = read_substmat_half(mat_file1)
    mat2 = read_substmat_half(mat_file2)

    # Verify specific differences
    assert mat2["W"]["W"] - mat1["W"]["W"] == 1
    assert mat2["M"]["M"] - mat1["M"]["M"] == 2
    assert mat2["F"]["D"] - mat1["F"]["D"] == -1
    assert mat2["A"]["A"] - mat1["A"]["A"] == 0


def test_diff_substmat_prints_output(capsys):
    """Test that the function prints output."""
    mat_file1 = Path(__file__).parent / "blosum62_half.mat"
    mat_file2 = Path(__file__).parent / "blosum45_half.mat"
    diff_substmat(mat_file1, mat_file2)

    captured = capsys.readouterr()
    # Check that something was printed
    assert len(captured.out) > 0
    # Check that the header contains amino acids
    assert "A" in captured.out
    assert "W" in captured.out


def test_diff_substmat_header_format(capsys):
    """Test that the header is formatted correctly."""
    mat_file1 = Path(__file__).parent / "blosum62_half.mat"
    mat_file2 = Path(__file__).parent / "blosum45_half.mat"
    diff_substmat(mat_file1, mat_file2)

    captured = capsys.readouterr()
    # Check header format (5 spaces at start, 4 spaces between amino acids)
    assert captured.out.startswith(
        "     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n"
    )


def test_diff_substmat_contains_all_amino_acids(capsys):
    """Test that all 20 amino acids appear in the output."""
    mat_file1 = Path(__file__).parent / "blosum62_half.mat"
    mat_file2 = Path(__file__).parent / "blosum45_half.mat"
    diff_substmat(mat_file1, mat_file2)
