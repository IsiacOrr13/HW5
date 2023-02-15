# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np


def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    x = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    G = x.align(hs_seq, gg_seq)
    M = x.align(hs_seq, mm_seq)
    B = x.align(hs_seq, br_seq)
    T = x.align(hs_seq, tt_seq)

    score_dict = dict()
    score_dict['Gallus gallus'] = G[0]
    score_dict['Mus musculus'] = M[0]
    score_dict['Balaeniceps_rex'] = B[0]
    score_dict['tursiops_truncatus'] = T[0]

    ordered_list_name = []
    ordered_list_score = []

    for i in range(4):
        key = max(score_dict, key=score_dict.get)
        score = score_dict[key]
        ordered_list_name.append(key)
        ordered_list_score.append(score)
        del score_dict[key]

    print('Species most similar to Homo Sapiens in Descending order: {q}'.format(q = ordered_list_name))
    print('Sequence Alignment Scores compared to Homo Sapiens in Descending order: {q}'.format(q = ordered_list_score))

    

if __name__ == "__main__":
    main()
