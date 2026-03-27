from create_substmat import create_substmat
from pathlib import Path
import pytest


def test_create_substmat_baba_output(capsys):
    """Test that baba.fasta produces the expected output."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    expected_output = (
        "     A    B    C\n" "A    1\n" "B   -2    3\n" "C   -2    1    3\n"
    )
    assert captured.out == expected_output


def test_create_substmat_header_format(capsys):
    """Test that the header is formatted correctly."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    # Check header format (5 spaces at start, 4 spaces between amino acids)
    assert captured.out.startswith("     A    B    C\n")


def test_create_substmat_diagonal_values(capsys):
    """Test that diagonal values are positive (self-matches)."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")

    # A-A should be 1
    assert "A    1" in lines[1]
    # B-B should be 3
    assert "B   -2    3" in lines[2]
    # C-C should be 3
    assert "C   -2    1    3" in lines[3]


def test_create_substmat_triangular_structure(capsys):
    """Test that output has triangular structure."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")

    # Skip header line
    data_lines = lines[1:]

    # First data line should have 1 value (A-A)
    assert data_lines[0].count("1") == 1
    # Second data line should have 2 values (B-A, B-B)
    assert len(data_lines[1].split()) == 3  # B + 2 values
    # Third data line should have 3 values (C-A, C-B, C-C)
    assert len(data_lines[2].split()) == 4  # C + 3 values


def test_create_substmat_invalid_amino_acid(capsys):
    """Test that invalid amino acids trigger error message."""
    mat_file = Path(__file__).parent / "baba.fasta"
    # Using default alphabet which doesn't include A, B, C
    create_substmat(mat_file)

    captured = capsys.readouterr()
    assert (
        "The alignment contains amino acids that are not in the provided alphabet"
        in captured.out
    )


def test_create_substmat_custom_alphabet(capsys):
    """Test with custom alphabet."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")

    # Should have header + 3 data lines
    assert len(lines) == 4
    # Header should only contain A, B, C
    header = lines[0].strip().split()
    assert set(header) == {"A", "B", "C"}


def test_create_substmat_output_not_empty(capsys):
    """Test that function produces output."""
    mat_file = Path(__file__).parent / "baba.fasta"
    create_substmat(mat_file, "ABC")

    captured = capsys.readouterr()
    assert len(captured.out) > 0
    assert captured.out.strip() != ""
