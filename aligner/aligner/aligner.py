from functools import partial, cache
import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from aligner import matrix_manipulation
from aligner.utilities import FastaEntry

BLOSUM62 = {
    'A': {'A':  4, 'C':  0, 'D': -2, 'E': -1, 'F': -2, 'G':  0, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S':  1, 'T':  0, 'V':  0, 'W': -3, 'Y': -2},
    'C': {'A':  0, 'C':  9, 'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1, 'M': -1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
    'D': {'A': -2, 'C': -3, 'D':  6, 'E':  2, 'F': -3, 'G': -1, 'H': -1, 'I': -3, 'K': -1, 'L': -4, 'M': -3, 'N':  1, 'P': -1, 'Q':  0, 'R': -2, 'S':  0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'E': {'A': -1, 'C': -4, 'D':  2, 'E':  5, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  1, 'L': -3, 'M': -2, 'N':  0, 'P': -1, 'Q':  2, 'R':  0, 'S':  0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F':  6, 'G': -3, 'H': -1, 'I':  0, 'K': -3, 'L':  0, 'M':  0, 'N': -3, 'P': -4, 'Q': -3, 'R': -3, 'S': -2, 'T': -2, 'V': -1, 'W':  1, 'Y':  3},
    'G': {'A':  0, 'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G':  6, 'H': -2, 'I': -4, 'K': -2, 'L': -4, 'M': -3, 'N':  0, 'P': -2, 'Q': -2, 'R': -2, 'S':  0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
    'H': {'A': -2, 'C': -3, 'D': -1, 'E':  0, 'F': -1, 'G': -2, 'H':  8, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N':  1, 'P': -2, 'Q':  0, 'R':  0, 'S': -1, 'T': -2, 'V': -3, 'W': -2, 'Y':  2},
    'I': {'A': -1, 'C': -1, 'D': -3, 'E': -3, 'F':  0, 'G': -4, 'H': -3, 'I':  4, 'K': -3, 'L':  2, 'M':  1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -2, 'T': -1, 'V':  3, 'W': -3, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D': -1, 'E':  1, 'F': -3, 'G': -2, 'H': -1, 'I': -3, 'K':  5, 'L': -2, 'M': -1, 'N':  0, 'P': -1, 'Q':  1, 'R':  2, 'S':  0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'L': {'A': -1, 'C': -1, 'D': -4, 'E': -3, 'F':  0, 'G': -4, 'H': -3, 'I':  2, 'K': -2, 'L':  4, 'M':  2, 'N': -3, 'P': -3, 'Q': -2, 'R': -2, 'S': -2, 'T': -1, 'V':  1, 'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F':  0, 'G': -3, 'H': -2, 'I':  1, 'K': -1, 'L':  2, 'M':  5, 'N': -2, 'P': -2, 'Q':  0, 'R': -1, 'S': -1, 'T': -1, 'V':  1, 'W': -1, 'Y': -1},
    'N': {'A': -2, 'C': -3, 'D':  1, 'E':  0, 'F': -3, 'G':  0, 'H':  1, 'I': -3, 'K':  0, 'L': -3, 'M': -2, 'N':  6, 'P': -2, 'Q':  0, 'R':  0, 'S':  1, 'T':  0, 'V': -3, 'W': -4, 'Y': -2},
    'P': {'A': -1, 'C': -3, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N': -2, 'P':  7, 'Q': -1, 'R': -2, 'S': -1, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
    'Q': {'A': -1, 'C': -3, 'D':  0, 'E':  2, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  1, 'L': -2, 'M':  0, 'N':  0, 'P': -1, 'Q':  5, 'R':  1, 'S':  0, 'T': -1, 'V': -2, 'W': -2, 'Y': -1},
    'R': {'A': -1, 'C': -3, 'D': -2, 'E':  0, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  2, 'L': -2, 'M': -1, 'N':  0, 'P': -2, 'Q':  1, 'R':  5, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
    'S': {'A':  1, 'C': -1, 'D':  0, 'E':  0, 'F': -2, 'G':  0, 'H': -1, 'I':  2, 'K':  0, 'L': -2, 'M': -1, 'N':  1, 'P': -1, 'Q':  0, 'R': -1, 'S':  4, 'T':  1, 'V': -2, 'W': -3, 'Y': -2},
    'T': {'A':  0, 'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N':  0, 'P': -1, 'Q': -1, 'R': -1, 'S':  1, 'T':  5, 'V':  0, 'W': -2, 'Y': -2},
    'V': {'A':  0, 'C': -1, 'D': -3, 'E': -2, 'F': -1, 'G': -3, 'H': -3, 'I':  3, 'K': -2, 'L':  1, 'M':  1, 'N': -3, 'P': -2, 'Q': -2, 'R': -3, 'S': -2, 'T':  0, 'V':  4, 'W': -3, 'Y': -1},
    'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F':  1, 'G': -2, 'H': -2, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y':  2},
    'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F':  3, 'G': -3, 'H':  2, 'I': -1, 'K': -2, 'L': -1, 'M': -1, 'N': -2, 'P': -3, 'Q': -1, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W':  2, 'Y':  7},
}

@cache
def blosum_62(aa1: str, aa2: str) -> int:
    """
    Looks up the score for two amino acids in the BLOSUM62 lookup table.

    If the amino acid code is not present in the lookup table, the lookup will
    raise a KeyError.

    Parameters
    ----------
    aa1 : str
    aa2 : str
        One-letter amino acid codes to be looked up.

    Returns
    -------
    int
        BLOSUM62 score.
    """
    return BLOSUM62[aa1.upper()][aa2.upper()]

def blosum_62_windowed(seq1: list[str], seq2: list[str]) -> int:
    """
    Returns the sum of BLOSUM62 scores for two lists of amino acids.

    If the amino acid code is not present in the lookup table, the lookup will
    raise a KeyError.

    Parameters
    ----------
    seq1 : [str]
    seq2 : [str]
        Lists of one-letter amino acid codes to be scored.

    Returns
    -------
    int
        BLOSUM62 score.
    """
    sum = 0
    for a, b in zip(seq1, seq2):
        sum += blosum_62(a, b)

    return sum

def hamming_distance(seq1: list[str], seq2: list[str]) -> int:
    """
    Returns the hamming distance between two sequences.

    Hamming distance is the number of unequal pairwise items in a sequence. In
    the example below, the Hamming distance is 3:
        AKLSACW
        |  | ||
        ALKSFCW

    The higher the Hamming distance, the lower the similarity between the
    sequences.

    Parameters
    ----------
    seq1 : [str]
    seq2 : [str]
        Sequences to be scored.

    Returns
    -------
    int
        Hamming distance between the two sequences.
    """
    dist = 0
    for a, b in zip(seq1, seq2):
        if a != b: dist += 1
        # Branchless variant:
        #   dist += int(a != b)

    return dist

def hamming_distance_inverse(seq1: list[str], seq2: list[str]) -> int:
    """
    Returns the inverse of the Hamming distance between two sequences. The
    Hamming distance is calculated normally, and is then subtracted from the
    length of the longest input sequence. This ensures that a fair score is
    returned even when the two sequences are not of equal length, making the
    function robust for edge cases.

    The inverse Hamming distance returns zero in the case that the sequences
    have no matching items (or if the shorter sequence has no matches in at 
    least the section of the shorter length in the longer sequence, in case the
    lengths are unequal).

    The higher the inverse Hamming distance, the greater the similarity. The
    maximum score returned is the length of the (shortest) input sequence.

    Parameters
    ----------
    seq1 : [str]
    seq2 : [str]
        Sequences to be scored.

    Returns
    -------
    int
        Inverse Hamming distance between the two sequences.
    """
    max = min(len(seq1), len(seq2))
    dist = hamming_distance(seq1, seq2)
    return max - dist

# From the Kyte & Doolittle article's C program:
# ```c
# char code[] "RKDBNSEHZQTGXAPVYCMILWF";
# float factor [] {0.0,0.6,1.0,1.0,1.0,3.6,1.0,1.3,1.0,1.0,3.8,4.1,4.1,6.3,
#     2.9,8.7,3.2,7.0,6.4,9.0,8.2,3.6,7.2};
# ```
#
# Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
# hydropathic character of a protein. Journal of molecular biology, 
# 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
KYTE_DOOLITTLE_FACTORS = {'R': 0.0, 'K': 0.6, 'D': 1.0, 'B': 1.0, 'N': 1.0,
                          'S': 3.6, 'E': 1.0, 'H': 1.3, 'Z': 1.0, 'Q': 1.0,
                          'T': 3.8, 'G': 4.1, 'X': 4.1, 'A': 6.3, 'P': 2.9,
                          'V': 8.7, 'Y': 3.2, 'C': 7.0, 'M': 6.4, 'I': 9.0,
                          'L': 8.2, 'W': 3.6, 'F': 7.2}

def kyte_doolittle_hydropathy_factor(aa: str) -> float:
    """
    Looks up the Kyte-Doolittle hydropathy factor for a given one-letter amino 
    acid code. The factor is returned as a float. 

    If the amino acid code is not present in the lookup table, the lookup will
    raise a KeyError.

    Parameters
    ----------
    aa : str
        A one-letter amino acid code.

    Returns
    -------
    float
        Kyte-Doolittle hydropathy factor for the provided amino acid.
    """
    return KYTE_DOOLITTLE_FACTORS[aa.upper()]

def hydrophobicity_score(aa1: str, aa2: str) -> int:
    """
    The higher, the more different the hydropathy between the input amino acids
    is. A score of zero corresponds to equal hydrophobicity according to the
    Kyle-Doolittle table. A high score suggests very dissimilar hydropathy.

    Parameters
    ----------
    aa1 : str
    aa2 : str
        One-letter amino acid codes to be compared.

    Returns
    -------
    int
        Hydropathy score calculated as |hydropathy_delta| * 10.
    """

    kd_aa1 = kyte_doolittle_hydropathy_factor(aa1)
    kd_aa2 = kyte_doolittle_hydropathy_factor(aa2)
    delta = abs(kd_aa1 - kd_aa2)

    return int(delta * 10)

def windowed(l: list, window_size, step=1) -> list[list]:
    """
    Returns a list of smaller lists with a equal to or smaller than
    `window_size`, each showing a window into `l` shifted by `step`. Once the
    number of remaining items in `l` is lower than `window_size`, the lists
    will decrease in length.

    The following example 
    `windowed([1, 2, 3, 4], 3)` 
    would return 
    `[[1, 2, 3], [2, 3, 4], [3, 4], [4]]`.

    Parameters
    ----------
    l : list
    window_size : int
    step : int, default 1

    Returns
    -------
    list[list]
        List containing the windows into `l`.
    """
    ret = []
    for i in range(len(l)):
        ret.append(l[i:i + window_size])

    return ret

def walk_up_down(M, pos: tuple[int, int]) -> tuple[
        tuple[int, int], tuple[int, int], int]:
    """
    Returns the start point, end point, and length of a particular non-zero
    diagonal stretch of a matrix by walking up and down from a given position.

    Note that the length of the diagonal is also returned, even though it is
    redundant for that information is already deducible from the start and end
    points. It is returned as a matter of convenience.

    Parameters
    ----------
    M : 2D numpy array
    pos : (int, int)

    Returns
    -------
    ((int, int), (int, int), int)
        A tuple with three items: the start position of the diagonal in M, the
        end position of the diagonal, and the length of the diagonal.
    """

    x, y = pos
    # The matrix is turned into a matrix of bools with True for non-zero cells,
    # and False otherwise. This makes detecting the end of the diagonal fast.
    M = M != 0

    # This convenience function returns True if the provided x and y
    # coordinates are within bounds of matrix M.
    shape = M.shape
    check_bounds = lambda x, y: 0 <= x < shape[1] and 0 <= y < shape[0]

    # First, we walk up until the upper zero or the edge is encountered.
    cur_x, cur_y = x, y
    while check_bounds(cur_x - 1, cur_y - 1) and M[cur_y - 1, cur_x - 1]:
        cur_y -= 1
        cur_x -= 1

    # The start point has now been encountered, and is stored.
    start_point = cur_x, cur_y
    
    # Now, we walk down until the lower zero or the edge is encountered.
    while check_bounds(cur_x + 1, cur_y + 1) and M[cur_y + 1, cur_x + 1]:
        cur_y += 1
        cur_x += 1

    # The end point has now been encountered, and is stored.
    end_point = cur_x, cur_y
    length = end_point[0] - start_point[0] + 1

    return start_point, end_point, length

class Aligner:
    """
    Aligner to study sequence homology between two sequences using a number of
    different alignment and filtering methods applied to a 2-dimensional matrix
    representation of the alignment.

    Attributes
    ----------
    s : FastaEntry
    t : FastaEntry
    window_size : int, default 5
    threshold : int, default 3
    alignment_matrix : 2D numpy array
        Matrix storing the alignment of `s` and `t`.
    filtering_matrix : 2D numpy array
        Matrix storing the filtered values derived from the alignment matrix.
    """

    def __init__(self, s: FastaEntry, t: FastaEntry, window_size=5, threshold=3):
        self.s = s
        self.t = t
        self.window_size = window_size
        self.threshold = threshold
        self._alignment_matrix = self._zeros()
        self._filtering_matrix = self.alignment_matrix
        self.alignment_matrix = self._alignment_matrix
        self.filtering_matrix = self._alignment_matrix

    @property
    def filtering_matrix(self):
        return self._filtering_matrix

    @filtering_matrix.setter
    def filtering_matrix(self, new_matrix):
        if self.filtering_matrix.shape == new_matrix.shape:
            self._filtering_matrix = new_matrix
        else:
            # new_matrix is not the same shape as self.filtering_matrix
            raise ValueError(
                "cannot set new filtering_matrix, because shapes do not match"
            )

    @property
    def alignment_matrix(self):
        return self._alignment_matrix

    @alignment_matrix.setter
    def alignment_matrix(self, new_matrix):
        if self.alignment_matrix.shape == new_matrix.shape:
            self._alignment_matrix = new_matrix
        else:
            # new_matrix is not the same shape as self.filtering_matrix
            raise ValueError(
                "cannot set new alignment_matrix, because shapes do not match"
            )

    def _zeros(self, dtype='int16'):
        """
        Convenience function for internal use which returns a 2D numpy array of
        the size of the alignment between `s` and `t` containing zeros of a
        `dtype`.

        Parameters
        ----------
        dtype, default 'int8'
            Datatype to initialize the matrix with.
        
        Returns
        -------
        np.array
        """
        return np.zeros((len(self.s.seq), len(self.t.seq)), dtype=dtype)

    def construct_alignment_matrix(self, scoring_func=lambda a, b: a == b): 
        """
        Constructs the single-point alignment matrix between `s` and `t`
        according to some `scoring_func` which by default places a 1 for any
        matches between the sequences, and otherwise 0.

        The resulting alignment matrix is returned, and the `alignment_matrix`
        and `filtering_matrix` attributes are set to the matrix.

        Parameters
        ----------
        scoring_func
            Function used to score alignment. Can return a number or a bool.
        
        Returns
        -------
        np.array
        """
        M = self._zeros()

        for i, a in enumerate(self.s.seq):
            for j, b in enumerate(self.t.seq):
                M[i, j] = scoring_func(a, b)

        self.alignment_matrix = M
        self.filtering_matrix = M

        return M

    def construct_windowed_alignment_matrix(self, scoring_func=hamming_distance_inverse):
        """
        Constructs the windowed alignment matrix between `s` and `t` according to some
        `alignment_func` which by default returns the inverse Hamming distance
        between the window sequences. The window size is taken from the `window_size` attribute.

        The resulting alignment matrix is returned, and the `alignment_matrix`
        and `filtering_matrix` attributes are set to the matrix.

        Parameters
        ----------
        alignment_func
            Function used to score alignment. Can return a number or a bool.
        
        Returns
        -------
        np.array
        """
        M = self._zeros()

        seq2_windowed = windowed(self.t.seq, self.window_size)
        for i, a in enumerate(windowed(self.s.seq, self.window_size)):
            for j, b in enumerate(seq2_windowed):
                score = scoring_func(a, b) 
                if score >= self.threshold:
                    M[i, j] = score

        self.alignment_matrix = M
        self.filtering_matrix = M

        return M

    def plot_matrix(self, M, title: str, truncate_label_length: int=32):
        """
        Shows a plot of matrix `M` with axes labeled according to the `id`
        field in the FastaEntry `s` and `t`.

        Parameters
        ----------
        M : np.array
            A 2D matrix representing data related to or contained in the
            Aligner class to be plotted.
        title : str
            Title of the plot.
        truncate_label_length : int, default
            The maximum number of characters to be displayed in the axis labels
            of the plot
        """
        fig, ax_with_seq = plt.subplots()
        ax_with_seq.set_title(title)
        ax_with_seq.set_xlabel(self.t.id[:truncate_label_length])
        ax_with_seq.set_ylabel(self.s.id[:truncate_label_length])
        ax_with_seq.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax_with_seq.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax_with_seq.matshow(M, interpolation='nearest')
        ax_with_seq.format_coord = lambda x, y: f"""x={int(x + 0.5)} \
({self.t.seq[int(x + 0.5)]}), \
y={int(y + 0.5)} \
({self.s.seq[int(y + 0.5)]})"""

        plt.xticks(range(len(self.t.seq)), self.t.seq)
        plt.yticks(range(len(self.s.seq)), self.s.seq)

        plt.show()

    def plot_alignment_matrix(self):
        """
        Shows a plot of the filtering matrix.
        """
        self.plot_matrix(self.alignment_matrix, "Alignment matrix")

    def plot_filter_matrix(self):
        """
        Shows a plot of the filtering matrix.
        """
        self.plot_matrix(self.filtering_matrix, "Filtering matrix")

    def longest_substrings(self) -> list[
            tuple[str, str, tuple[int, int], tuple[int, int], int]]:
        """
        Returns a list of the matching longest substrings according to the
        alignment matrix and filtering matrix by searching for long diagonals
        by finding maximum-value cells.

        The positions of cells which have the maximum value of the
        filtering_matrix are stored. For each of these positions, the diagonal
        is walked up and down until the diagonal terminates, and the resulting
        start position, end position, and diagonal length are stored in a set.
        The longest substrings are then retrieved from `s` and `t`.

        Returns
        -------
        list[(str, str, (int, int), (int, int), int)]
            A sorted list of 5-item tuple containing the `s` and `t` substring,
            the start and end points, and the length.
        """
        # Determine the matrix indicese of the cells in the filtering matrix
        # which are equal to the maximum value in the matrix.
        maxes = self.filtering_matrix == self.filtering_matrix.max()
        # Get positions of these maximum-value points.
        ys, xs = np.asarray(maxes).nonzero()

        results = set()
        for y, x in zip(ys, xs):
            res = walk_up_down(self.alignment_matrix, (x, y))
            results.add(res)

        longest_substrings = [("".join(self.s.seq[start[1]:end[1]+1]), 
                               "".join(self.t.seq[start[0]:end[0]+1]), 
                               start, end, length) 
                              for start, end, length in results]
        # Sort according to substring length as specified in the 5th field.
        longest_substrings.sort(key=lambda entry: entry[4])

        return longest_substrings

    def _generic_pass(self, pass_function):
        """
        An generic pass function for internal use. Takes a pass function to
        apply to the filtering matrix, and sets the `filtering_matrix`
        attribute to the value returned by `pass_function`. That value is also
        returned by the function.

        Parameters
        ----------
        pass_function
            A function that can take a 2D matrix such as the filtering matrix
            and return some matrix of the same shape.

        Returns
        -------
        np.array
            The new filtering_matrix, as returned by `pass_function`.
        """
        a = pass_function(self.filtering_matrix)
        self.filtering_matrix = a
        return a

    def addition_pass(self):
        """
        Apply an addition pass over the filtering matrix.

        Within existing non-zero values, the values diagonally above and below
        every cell are added to that cells value. High values will emerge
        towards the center of diagonals, allowing for identification of long
        diagonals.

        Returns
        -------
        np.array
            The filtering matrix after application of an addition pass.
        """
        return self._generic_pass(matrix_manipulation.addition_pass)

    def shortening_pass(self):
        """
        Apply a shortening pass over the filtering matrix.

        All cells are multiplied with the values in their diagonally adjacent
        cells. If one of these adjacent cells contains 0, this cell will also
        be set to 0. Multiple of these passes will result in shortening of
        diagonals, and elimination of short diagonals.

        Returns
        -------
        np.array
            The filtering matrix after application of a shortening pass.
        """
        return self._generic_pass(matrix_manipulation.shortening_pass)

    def normalize_pass(self, to=1):
        """
        Apply a normalization pass over the filtering matrix.

        All non-zero values in the filtering matrix are normalized to a single
        value, `to` (default=1).

        Parameters
        ----------
        to : int, default 1
            Value to which any non-zero value is normalized

        Returns
        -------
        np.array
            The filtering matrix after application of a normalization pass.
        """
        return self._generic_pass(matrix_manipulation.normalize_pass)

    def threshold_pass(self, threshold=None):
        """
        Apply a threshold pass over the filtering matrix.

        Any value below a threshold is set to zero. The default behavior is to
        use the `threshold` attribute of the class, unless it is specified when
        called.

        Parameters
        ----------
        threshold : int, default None
            Threshold value under which cells are set to zero. If None
            (default), the class attribute of `threshold` will be used.

        Returns
        -------
        np.array
            The filtering matrix after application of a threshold pass.
        """
        if threshold == None: 
            thr = self.threshold
        else: 
            thr = threshold
        return self._generic_pass(
            partial(matrix_manipulation.threshold_pass, threshold=thr)
        )

    def squaring_pass(self):
        """
        Apply a squaring pass over the filtering matrix.

        The values in the matrix are squared.

        Returns
        -------
        np.array
            The filtering matrix after application of a squaring pass.
        """
        return self._generic_pass(matrix_manipulation.squaring_pass)

    def log_pass(self):
        """
        Apply a logarithm pass over the filtering matrix.

        The logarithm of the values in the matrix is taken.

        Returns
        -------
        np.array
            The filtering matrix after application of a logarithm pass.
        """
        return self._generic_pass(matrix_manipulation.log_pass)

    def lifting_pass(self):
        """
        Apply a lifting pass over the filtering matrix.

        The values in the matrix are lifted, by adding the absolute value of
        the minimum value in the matrix to every cell. If the matrix contains a
        negative minimum value, all values in the matrix will be greater or
        equal to zero after a lifting pass.

        Returns
        -------
        np.array
            The filtering matrix after application of a lifting pass.
        """
        return self._generic_pass(matrix_manipulation.lifting_pass)

    def long_substrings_report(self):
        """
        Find long substrings by 

        The values in the matrix are lifted, by adding the absolute value of
        the minimum value in the matrix to every cell. If the matrix contains a
        negative minimum value, all values in the matrix will be greater or
        equal to zero after a lifting pass.

        Returns
        -------
        np.array
            The filtering matrix after application of a lifting pass.
        """
        print("Comparing")
        print(f"\ts:  {len(self.s.seq)}\t{self.s.id}")
        print(f"\tt:  {len(self.t.seq)}\t{self.t.id}")

        # Time this endeavor.
        t_start = time.perf_counter()

        # From the filtering matrix, strong diagonals can be resolved, which
        # point to the longest substrings which might show interesting
        # homologies.
        ls = self.longest_substrings()

        t_end = time.perf_counter()
        total_time = round(t_end - t_start, 3)
        print(f"Done, completed in {total_time} seconds")

        # Print a nice header for the of substrings and their positions.
        print(f"\t{'start':<10}\t{'end':<10}\tlen ")
        print(f"\t{'—' * 10}\t{'—' * 10}\t{'—' * 4}") 
        # Print the top 10 longest substrings.
        # (Using this dreadful-looking f-string and list comprehension-join
        # apperatus.)
        # This diff function returns a string of ' ' and '|' signifying a
        # mismatch or a match for each position in the substrings,
        # respectively. The [...][expr] bit is clever but not very readable 
        # i guess.
        diff = lambda s_sub, t_sub: "".join([
            [" ", "|"][s == t] for s, t in zip(list(s_sub), list(t_sub))])
        print("\n".join(
            [f"""\t{str(start):<10}\t{str(end):<10}\t{length:>4}
s:  {s_substring}
    {diff(s_substring, t_substring)}
t:  {t_substring}"""
             for s_substring, t_substring, start, end, length in ls[::-1][:10]]
        ))
