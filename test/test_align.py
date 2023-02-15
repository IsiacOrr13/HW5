# Importing Dependencies
import pytest
from align import read_fasta, NeedlemanWunsch

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    x = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    x.align(seq1, seq2)
    align_mat = x._align_matrix
    align_vec = x._align_matrix_vector
    assert align_mat[2] == [-2.0, 0.0, 0.0, 1.0, 0.0]
    assert align_vec[1] == [0, 1, 2, 2, 2]
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    x = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    y = x.align(seq3, seq4)
    assert y[1] == 'MAVHQLIRRP'
    assert y[2] == 'M---QLIRHP'





