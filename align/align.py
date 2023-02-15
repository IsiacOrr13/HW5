# Importing Dependencies
import numpy as np
from typing import Tuple


# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._align_matrix_vector = None
        self._gapA = 0
        self._gapB = 0

        self.M = 1
        self.m = -1
        self.g_ex = gap_extend
        self.g_open = gap_open

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm

        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA

        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """

        self._align_matrix = [[0] * (len(seqA)+1) for i in range((len(seqB)+1))]
        self._align_matrix_vector = [[0] * (len(seqA)+1) for i in range((len(seqB)+1))]

        for h in range(len(seqA)+1):
            self._align_matrix[0][h] = self.g_ex * h
        for k in range(len(seqB)+1):
            self._align_matrix[k][0] = self.g_ex * k

        for i in range(len(seqB)):
            i += 1
            for j in range(len(seqA)):
                j += 1

                up = self._align_matrix[i - 1][j] + self.g_ex

                back = self._align_matrix[i][j - 1] + self.g_ex

                if seqB[i - 1] == seqA[j - 1]:
                    diag = self._align_matrix[i - 1][j - 1] + self.M
                else:
                    diag = self._align_matrix[i - 1][j - 1] + self.m

                self._align_matrix[i][j] = max(up, diag, back)

                vecs = [up, diag, back]
                vec = np.argmax(vecs)
                self._align_matrix_vector[i][j] = vec

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        return self._backtrace(seqA, seqB)

    def _backtrace(self, seqA, seqB) -> Tuple[float, str, str]:
        """
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        i = -1
        j = -1
        while abs(j) <= len(seqA) and abs(i) <= len(seqB):
            if self._align_matrix_vector[i][j] == 1:
                self.seqA_align = seqA[j] + self.seqA_align
                self.seqB_align = seqB[i] + self.seqB_align
                self.alignment_score += (self.sub_dict[(seqA[j], seqB[i])])
                if self._gapA != 0:
                    self.alignment_score += (self.g_open + (self.g_ex * self._gapA))
                if self._gapB != 0:
                    self.alignment_score += (self.g_open + (self.g_ex * self._gapB))
                self._gapA = 0
                self._gapB = 0
                i -= 1
                j -= 1

            elif self._align_matrix_vector[i][j] == 0:
                self.seqA_align = "-" + self.seqA_align
                self.seqB_align = seqB[i] + self.seqB_align
                i -= 1
                self._gapA += 1
            else:
                self.seqA_align = seqA[j] + self.seqA_align
                self.seqB_align = "-" + self.seqB_align
                j -= 1
                self._gapB += 1

        return (self.alignment_score, self.seqA_align, self.seqB_align)

def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header














