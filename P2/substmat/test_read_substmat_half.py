from pathlib import Path
from read_substmat_half import read_substmat_half


def test_read_substmat_minimal(tmp_path):
    """
    Test read_substmat using a substitution matrix file.
    """

    mat_file = Path(__file__).parent / "blosum62_half.mat"

    # Read the matrix using the function
    smat = read_substmat_half(str(mat_file))

    #  Checks equivalent to the original doctest
    assert smat["W"]["W"] == 11
    assert smat["W"]["A"] == -3
    assert smat["A"]["W"] == -3
    assert smat["A"]["A"] == 4

    # Check dimensions
    assert len(smat) == 20
    assert len(smat["W"]) == 20
    assert len(smat["A"]) == 20

    #  Check presence of keys
    assert "W" in smat["A"]
    assert "A" in smat["W"]
