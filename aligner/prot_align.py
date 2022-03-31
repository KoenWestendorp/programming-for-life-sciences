from collections import namedtuple
import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import signal

from main import addition_pass, normalize_pass, shortening_pass
import utilities

FastaEntry = namedtuple('FastaEntry', 'id seq')

def construct_alignment_matrix(seq1, seq2):
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            if a == b: M[i, j] = 1

    return M

def construct_amino_acid_matrix(seq1, seq2):
    assert False
    M = np.array((len(seq1), len(seq2)), dtype=tuple[str])

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            M[i, j] = str()

    return M

def hamming_distance(seq1: list[str], seq2: list[str]) -> int:
    dist = 0
    for a, b in zip(seq1, seq2):
        if a != b: dist += 1

    return dist

def hamming_distance_inverse(seq1: list[str], seq2: list[str]) -> int:
    max = len(seq1)
    dist = hamming_distance(seq1, seq2)
    return max - dist

def windowed(l: list, window_size) -> list:
    ret = []
    for i in range(len(l) - window_size):
        ret.append(l[i:i + window_size])

    return ret

def construct_windowed_alignment_matrix(seq1, seq2, n=5, scoring_func=hamming_distance_inverse, scoring_threshold=1):
    # TODO: make this function take a parameter for distance function, with
    # Hamming distance function as default.
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    for i, a in enumerate(windowed(seq1, n)):
        for j, b in enumerate(windowed(seq2, n)):
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
        x < shape[0] and y < shape[1] and x >= 0 and y >= 0

    # First, we walk up until the upper zero or the edge is encountered.
    yt, xt = int(y), int(x)
    while check_bounds(M.shape, xt - 1, yt - 1) and M[yt - 1, xt - 1]:
        yt -= 1
        xt -= 1

    start_point = xt, yt
    
    while check_bounds(M.shape, xt + 1, yt + 1) and M[yt + 1, xt + 1]:
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

def compare(s: FastaEntry, t: FastaEntry):
    print("Comparing")
    print(f"\t{len(s.seq)}\t{s.id}")
    print(f"\t{len(t.seq)}\t{t.id}")

    t_start = time.perf_counter()

    M = construct_windowed_alignment_matrix(s.seq, t.seq, n=5)
    N = M.copy()

    print(f"Matrix constructed ({round(time.perf_counter() - t_start, 3)} seconds)")

    M = addition_pass(M)
    A = M.copy()

    M = addition_pass(M)
    M = normalize_pass(M)
    M = addition_pass(M)
    M = shortening_pass(M)

    print(f"Passes completed ({round(time.perf_counter() - t_start, 3)} seconds)")

    ls = longest_substrings(M, N, s)

    t_end = time.perf_counter()

    print(f"Done, completed in {round(t_end - t_start, 3)} seconds")

    print(f"\tlen \t{'start':<10}\t{'end':<10}\tsequence")
    print("\n".join([f"\t{end[0] - start[0]:>4}\t{str(start):<10}\t{str(end):<10}\t{seq}" for seq, start, end in ls[::-1][:10]]))

    plot_matrix(A, s, t)

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

        compare(entries[0], entries[1])
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

