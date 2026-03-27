from pathlib import Path
from read_substmat_full import read_substmat_full


def test_read_substmat_full():

    mat_file = Path(__file__).parent / "blosum62_full.mat"
    
    smat = read_substmat_full(str(mat_file))

    assert smat["W"]["W"] == 11
    assert smat["W"]["A"] == -3
    assert smat["A"]["W"] == -3
    assert len(smat) == 20
    assert len(smat["L"]) == 20
