from collections import namedtuple
from itertools import islice

import matplotlib.pyplot as plt
import numpy as np

from main import addition_pass, normalize_pass, shortening_pass
import utilities

FastaEntry = namedtuple('FastaEntry', 'id seq')

def construct_alignment_matrix(seq1, seq2):
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            if a == b: M[i, j] = 1

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

def construct_windowed_alignment_matrix(seq1, seq2, n=5, scoring_func=hamming_distance_inverse):
    # TODO: make this function take a parameter for distance function, with
    # Hamming distance function as default.
    M = np.zeros((len(seq1), len(seq2)), dtype='int8')

    for i, a in enumerate(windowed(seq1, n)):
        for j, b in enumerate(windowed(seq2, n)):
            M[i, j] = scoring_func(a, b)

    return M

def windowed(l: list, window_size) -> list:
    ret = []
    for i in range(len(l) - window_size):
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

def main(s, t):
    fig = plt.figure()

    M = construct_windowed_alignment_matrix(s, t, n=5)
    N = M.copy()

    ax1 = fig.add_subplot(131, title="Windowed alignment matrix (n=5)")
    ax1.imshow(M, interpolation='nearest')

    M = addition_pass(M)

    ax2 = fig.add_subplot(132, title="After 1 addition pass")
    ax2.imshow(M, interpolation='nearest')

    fig_with_seq = plt.figure()
    ax_with_seq = fig_with_seq.add_subplot(111, 
                                           xlabel="SARS-CoV-2 ORF3a", 
                                           xticks=range(len(t)), 
                                           ylabel="SARS-CoV-1 M protein", 
                                           yticks=range(len(s)))
    ax_with_seq.set_xticklabels(t, fontdict={'fontsize': 9})
    ax_with_seq.set_yticklabels(s, fontdict={'fontsize': 9})
    ax_with_seq.imshow(M, interpolation='bilinear')

    M = addition_pass(M)
    M = addition_pass(M)
    M = normalize_pass(M)
    M = addition_pass(M)
    M = addition_pass(M)
    M = shortening_pass(M)
    M = normalize_pass(M)

    ax3 = fig.add_subplot(133, title="After 3 addition passes, \nnormalization, \n2 more addition passes")
    ax3.imshow(M, interpolation='nearest')

    ys, xs = np.asarray(M == M.max()).nonzero()
    #print(xs, ys)
    results = set()
    for y, x in zip(ys, xs):
        res = walk_up_down(N, (x, y))
        start_point, end_point, diagonal_len = res
        results.add(res)

        #print(f"diagonal through {x, y}: {start_point}, {end_point}, {diagonal_len}")
        #print(results)

    longest_substrings = [("".join(s[start_point[1]:end_point[1]+1]), start_point, end_point) 
                          for start_point, end_point, _ 
                          in results]

    print(longest_substrings)

    plt.show()
    #plt.savefig("figure.svg")

if __name__ == "__main__":
    s1, *_ = utilities.read_fasta_file("P59596.fasta")
    sars_cov_1_m_protein = FastaEntry(s1[0], s1[1])
    s2, *_ = utilities.read_fasta_file("P0DTC5.fasta")
    sars_cov_2_m_protein = FastaEntry(s2[0], s2[1])
    s3, *_ = utilities.read_fasta_file("P0DTC3.fasta")
    sars_cov_2_orf3a = FastaEntry(s3[0], s3[1])

    s = sars_cov_1_m_protein.seq 
    t = sars_cov_2_orf3a.seq

    #s = "12345abc67890"
    #t = ".',<>abc_/=[]"
    #s = "eoi12345abc67890trl"
    #t = "eoi.',<>abc_/=[]trl"

    main(s, t)

