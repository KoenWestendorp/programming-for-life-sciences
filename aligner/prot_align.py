from collections import namedtuple
import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import signal

from main import addition_pass, normalize_pass, shortening_pass
import utilities

FastaEntry = namedtuple('FastaEntry', 'id seq')

def construct_alignment_matrix(seq1, seq2, scoring_func=lambda a, b: a == b):
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            if scoring_func: M[i, j] = 1

    return M

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

def construct_windowed_alignment_matrix(seq1, seq2, n=5, scoring_func=hamming_distance_inverse, scoring_threshold=1):
    # TODO: make this function take a parameter for distance function, with
    # Hamming distance function as default.
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    seq2_windowed = windowed(seq2, n)
    for i, a in enumerate(windowed(seq1, n)):
        for j, b in enumerate(seq2_windowed):
            score = scoring_func(a, b) 
            if score >= scoring_threshold:
                M[i, j] = score

    return M

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

def plot_matrix(M, s, t):
    fig_with_seq, ax_with_seq = plt.subplots()
    ax_with_seq.set_xlabel(t.id)
    ax_with_seq.set_ylabel(s.id)
    #ax_with_seq.set_xticks(range(len(t.seq)), minor=True)
    #ax_with_seq.set_yticks(range(len(s.seq)), minor=True)
    ax_with_seq.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax_with_seq.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax_with_seq.matshow(M, interpolation='nearest')
    ax_with_seq.format_coord = lambda x, y: f"x={int(x+0.5)} ({t.seq[int(x+0.5)]}), y={int(y+0.5)} ({s.seq[int(y+0.5)]})"

    plt.xticks(range(len(t.seq)), t.seq)
    plt.yticks(range(len(s.seq)), s.seq)

    plt.show()

def longest_substrings(M, N, s):
    ys, xs = np.asarray(M == M.max()).nonzero()
    results = set()
    for y, x in zip(ys, xs):
        res = walk_up_down(N, (x, y))
        #start_point, end_point, diagonal_len = res
        results.add(res)

    longest_substrings = [("".join(s.seq[start_point[1]:end_point[1]+1]), start_point, end_point) 
                          for start_point, end_point, _ 
                          in results]

    longest_substrings.sort(key=lambda e: e[2][0] - e[1][0] + 1)

    return longest_substrings

def compare(s: FastaEntry, t: FastaEntry, n=5, threshold=2):
    print("Comparing")
    print(f"\t{len(s.seq)}\t{s.id}")
    print(f"\t{len(t.seq)}\t{t.id}")

    t_start = time.perf_counter()

    M = construct_alignment_matrix(s.seq, t.seq, scoring_func=blosum_62)
    #M = construct_windowed_alignment_matrix(s.seq, t.seq, scoring_func=blosum_62_windowed)
    N = M.copy()

    print(f"Matrix constructed ({round(time.perf_counter() - t_start, 3)} seconds)")

    M = addition_pass(M)
    M = addition_pass(M)
    M = addition_pass(M)
    A = M.copy()

    M = addition_pass(M)
    M = normalize_pass(M)
    M = addition_pass(M)
    M = shortening_pass(M)
    M[M <= threshold] = 0
    #M += N

    print(f"Passes completed ({round(time.perf_counter() - t_start, 3)} seconds)")

    ls = longest_substrings(M, N, s)

    t_end = time.perf_counter()

    print(f"Done, completed in {round(t_end - t_start, 3)} seconds")

    print(f"\tlen \t{'start':<10}\t{'end':<10}\tsequence")
    print("\n".join([f"\t{end[0] - start[0]:>4}\t{str(start):<10}\t{str(end):<10}\t{seq}" for seq, start, end in ls[::-1][:10]]))


    plot_matrix(N, s, t)
    plot_matrix(A, s, t)
    return M

def main():
    import sys

    if len(sys.argv) == 3:
        fasta_files = sys.argv[1:3]

        # Assuming single-entry fasta files
        #entries = [utilities.read_fasta_file(file_path)[0] for file_path in fasta_files]
        entries = []
        for file_path in fasta_files:
            print(f"Reading '{file_path}'...")
            entry = utilities.read_fasta_file(file_path)[0]
            entries.append(FastaEntry(entry[0], entry[1]))

        #t_max = 3
        #fig, axs = plt.subplots(1, t_max)
        #for t in range(0, t_max):
        #    #for ni, n in enumerate(range(1, 4)):
        #    n = 2
        #    M = compare(entries[0], entries[1], n=n, threshold=t)
        #    axs[t].matshow(M, interpolation='nearest')

        #plt.show()

        compare(entries[0], entries[1], n=3, threshold=3)

    else:
        print("Please provide two paths to fasta files containing protein sequences you want to compare.")

if __name__ == "__main__":
    main()

    # s1, *_ = utilities.read_fasta_file("P59596.fasta")
    # sars_cov_1_m_protein = FastaEntry(s1[0], s1[1])
    # s2, *_ = utilities.read_fasta_file("P0DTC5.fasta")
    # sars_cov_2_m_protein = FastaEntry(s2[0], s2[1])
    # s3, *_ = utilities.read_fasta_file("P0DTC3.fasta")
    # sars_cov_2_orf3a = FastaEntry(s3[0], s3[1])

    # s = sars_cov_2_m_protein.seq 
    # t = sars_cov_2_orf3a.seq

    # compare(s, t)

