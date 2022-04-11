"""
Functions used for manipulating alignment matrices (numpy arrays).

The documentation in this file is sparse because the logic, I believe, is
rather trivial, or has already been described in detail in the documentation
for the methods of aligner.Aligner. I rather present a single source of truth
at the more user-facing side, than keeping track of two descriptions which
might get out of sync. Therefore, explanation of the operation of the functions
implemented here can be found in `aligner.py`, at the respective
aligner.Aligner methods.
"""

import numpy as np

def add_diagonals(N):
    h, w = N.shape
    R = N.copy()
    for y in range(h):
        for x in range(w):
            try:
                diag_above = N[y - 1, x - 1]
            except:
                # Out of bounds
                diag_above = 0

            try:
                diag_below = N[y + 1, x + 1]
            except:
                # Out of bounds
                diag_below = 0

            r = 0

            if diag_above > 0: r += diag_above
            if diag_below > 0: r += diag_below
            if N[y, x] == 0: r = 0

            R[y, x] = r

    return R

def shorten_diagonals(N):
    h, w = N.shape
    R = N.copy()
    for y in range(h):
        for x in range(w):
            try:
                diag_above = N[y - 1, x - 1]
            except:
                # Out of bounds
                diag_above = 0

            try:
                diag_below = N[y + 1, x + 1]
            except:
                # Out of bounds
                diag_below = 0

            r = N[y, x]

            r *= diag_above
            r *= diag_below

            R[y, x] = r

    return R

def addition_pass(M):
    M = add_diagonals(M)
    return M

def shortening_pass(M):
    M = shorten_diagonals(M)
    return M

def normalize_pass(M, to=1):
    M[M > 0] = to
    return M 

def threshold_pass(M, threshold: int=1):
    M[M <= threshold] = 0
    return M

def squaring_pass(M):
    return M ** 2

def log_pass(M):
    return np.log(M)

def lifting_pass(M):
    lowest = M.min()
    return M + abs(lowest)

def print_matrix_properties(M, print_matrix=True):
    if print_matrix: print(M)
    print("  max:", M.max())
    print("n_max:", np.count_nonzero(M == M.max()))
    print(" ones:", np.count_nonzero(M == 1))
    print("shape:", M.shape)
