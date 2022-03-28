# Koen Westendorp
# 2022-03-28

import numpy as np

def add_diagonals(N):
    # for y, row in enumerate(N):
    #     for x, cell in enumerate(row):
    #         if 

    h, w = N.shape
    R = N.copy()
    for y in range(h):
        for x in range(w):
            try:
                diag_above = N[y - 1, x - 1]
            except:
                diag_above = 0

            try:
                diag_below = N[y + 1, x + 1]
            except:
                diag_below = 0

            r = 0

            if diag_above > 0: r += diag_above
            if diag_below > 0: r += diag_below
            if N[y, x] == 0: r = 0

            R[y, x] = r

    return R

def pad(N):
    return np.pad(N, (1, 1), 'constant', constant_values=(0, 0))

def add_diagonals_pad(N):
    return pad(add_diagonals(N))

if __name__ == "__main__":
    s = list("actgatgcaag")
    t = list("ctgatgcata")
    # s = list("hactgu")
    # t = list("dactgo")

    M = np.zeros((len(s), len(t)))

    for i, a in enumerate(s):
        for j, b in enumerate(t):
            if a == b: M[i, j] = 1

    print("Comparison")
    print(M)

    for n in range(6):
        print("Pass", n)

        M = add_diagonals(M)

        print(M)
        print(M.shape)
