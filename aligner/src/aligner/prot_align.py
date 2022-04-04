from collections import namedtuple
from functools import partial
import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from main import addition_pass, normalize_pass, shortening_pass, threshold_pass, squaring_pass, lifting_pass
import utilities

FastaEntry = namedtuple('FastaEntry', 'id seq')

def blosum_62(a: str, b: str) -> int:
    matrix = {
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

    return matrix[a][b]

def blosum_62_windowed(seq1: list[str], seq2: list[str]) -> int:
    sum = 0
    for a, b in zip(seq1, seq2):
        sum += blosum_62(a, b)

    return sum

def hamming_distance(seq1: list[str], seq2: list[str]) -> int:
    dist = 0
    for a, b in zip(seq1, seq2):
        if a != b: dist += 1

    return dist

def hamming_distance_inverse(seq1: list[str], seq2: list[str]) -> int:
    max = min(len(seq1), len(seq2))
    dist = hamming_distance(seq1, seq2)
    return max - dist

def windowed(l: list, window_size) -> list:
    ret = []
    for i in range(len(l)):
        ret.append(l[i:i + window_size])

    return ret

def walk_up_down(M, pos: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int], int]:
    """
    Returns the start point, end point, and length of a particular non-zero
    diagonal stretch of a matrix by walking up and down from a given position.
    """

    x, y = pos
    M = M != 0

    check_bounds = lambda shape, x, y: \
        x < shape[1] and y < shape[0] and x >= 0 and y >= 0

    # First, we walk up until the upper zero or the edge is encountered.
    yt, xt = int(y), int(x)
    while check_bounds(M.shape, xt - 1, yt - 1) and M[yt - 1, xt - 1]:
        yt -= 1
        xt -= 1

    start_point = xt, yt
    
    while check_bounds(M.shape, xt + 1, yt + 1)\
        and M[yt + 1, xt + 1]:
        yt += 1
        xt += 1

    end_point = xt, yt
    length = end_point[0] - start_point[0] + 1

    return start_point, end_point, length

class Aligner:
    def __init__(self, s: FastaEntry, t: FastaEntry, window_size=5, threshold=3):
        self.s = s
        self.t = t
        self.window_size = window_size
        self.threshold = threshold
        self._alignment_matrix = self._zeros()
        self._filtering_matrix = self._alignment_matrix

    @property
    def filtering_matrix(self):
        return self._filtering_matrix

    @filtering_matrix.setter
    def filtering_matrix(self, new_matrix):
        if self.filtering_matrix.shape == new_matrix.shape:
            self._filtering_matrix = new_matrix
        else:
            # new_matrix is not the same shape as self.filtering_matrix
            # TODO: Find a more appropriate exception for this case.
            raise TypeError

    def _zeros(self, dtype='int8'):
        return np.zeros((len(self.s.seq), len(self.t.seq)), dtype=dtype)

    def construct_alignment_matrix(self, scoring_func=lambda a, b: a == b):
        M = self._zeros()

        for i, a in enumerate(self.s.seq):
            for j, b in enumerate(self.t.seq):
                M[i, j] = scoring_func(a, b)

        self.alignment_matrix = M

        return M

    def construct_windowed_alignment_matrix(self, scoring_func=hamming_distance_inverse):
        # TODO: make this function take a parameter for distance function, with
        # Hamming distance function as default.
        M = self._zeros()

        seq2_windowed = windowed(self.t.seq, self.window_size)
        for i, a in enumerate(windowed(self.s.seq, self.window_size)):
            for j, b in enumerate(seq2_windowed):
                score = scoring_func(a, b) 
                if score >= self.threshold:
                    M[i, j] = score

        self.alignment_matrix = M

        return M

    def plot_matrix(self, M):
        _, ax_with_seq = plt.subplots()
        ax_with_seq.set_xlabel(self.t.id)
        ax_with_seq.set_ylabel(self.s.id)
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
        self.plot_matrix(self.alignment_matrix)

    def plot_filter_matrix(self):
        self.plot_matrix(self.filtering_matrix)

    def longest_substrings(self):
        ys, xs = np.asarray(self.filtering_matrix 
                            == self.filtering_matrix.max()).nonzero()
        results = set()
        for y, x in zip(ys, xs):
            res = walk_up_down(self.alignment_matrix, (x, y))
            results.add(res)

        longest_substrings = [("".join(self.s.seq[start[1]:end[1]+1]), 
                               start, end) 
                              for start, end, _ in results]
        longest_substrings.sort(key=lambda e: e[2][0] - e[1][0] + 1)

        return longest_substrings

    def _generic_pass(self, pass_function):
        a = pass_function(self.filtering_matrix)
        self.filtering_matrix = a
        return a

    def addition_pass(self):
        return self._generic_pass(addition_pass)

    def shortening_pass(self):
        return self._generic_pass(shortening_pass)

    def normalize_pass(self):
        return self._generic_pass(normalize_pass)

    def threshold_pass(self):
        return self._generic_pass(
            partial(threshold_pass, threshold=self.threshold)
        )

    def squaring_pass(self):
        return self._generic_pass(squaring_pass)

    def lifting_pass(self):
        return self._generic_pass(lifting_pass)

    def compare(self):
        print("Comparing")
        print(f"\t{len(self.s.seq)}\t{self.s.id}")
        print(f"\t{len(self.t.seq)}\t{self.t.id}")

        t_start = time.perf_counter()

        self.construct_windowed_alignment_matrix(scoring_func=blosum_62_windowed)

        print(f"Matrix constructed ({round(time.perf_counter() - t_start, 3)} seconds)")

        self.lifting_pass()
        self.addition_pass()
        self.lifting_pass()
        self.addition_pass()

        print(f"Passes completed ({round(time.perf_counter() - t_start, 3)} seconds)")

        ls = self.longest_substrings()

        t_end = time.perf_counter()

        print(f"Done, completed in {round(t_end - t_start, 3)} seconds")

        print(f"\tlen \t{'start':<10}\t{'end':<10}\tsequence")
        print("\n".join([f"\t{end[0] - start[0]:>4}\t{str(start):<10}\t{str(end):<10}\t{seq}" for seq, start, end in ls[::-1][:10]]))

        self.plot_alignment_matrix()
        self.plot_filter_matrix()

        return self.filtering_matrix

def main():
    import sys

    if len(sys.argv) == 3:
        fasta_files = sys.argv[1:3]

        # Assuming single-entry fasta files
        entries = []
        for file_path in fasta_files:
            print(f"Reading '{file_path}'...")
            entry = utilities.read_fasta_file(file_path)[0]
            entries.append(FastaEntry(entry[0], entry[1]))

        al = Aligner(entries[0], entries[1], window_size=3, threshold=3)
        al.compare()

    else:
        print("Please provide two paths to fasta files containing protein sequences you want to compare.")

if __name__ == "__main__":
    main()
