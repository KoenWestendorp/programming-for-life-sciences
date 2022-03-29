# Koen Westendorp
# 2022-03-28

import numpy as np

def add_diagonals(N):
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

def shorten_diagonals(N):
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

            r = N[y, x]

            r *= diag_above
            r *= diag_below

            R[y, x] = r

    return R

def addition_pass(M):
    M = add_diagonals(M)

    # print(M)
    # print("  max:", M.max())
    # print("n_max:", np.count_nonzero(M == M.max()))
    # print(" ones:", np.count_nonzero(M == 1))
    # print("shape:", M.shape)

    return M

def shortening_pass(M):
    M = shorten_diagonals(M)

    return M

def normalize_pass(M):
    M[M > 0] = 1
    return M 

def walk_up_down(M, pos: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int], int]:
    """
    Returns the start point, end point, and length of a particular non-zero
    diagonal stretch of a matrix by walking up and down from a given position.
    """

    x, y = pos
    M = M != 0.

    # First, we walk up until the upper zero or the edge is encountered.
    yt, xt = int(y), int(x)
    while M[yt, xt]:
        yt -= 1
        xt -= 1

    # This overcounts the xt and yt values by one, because the condition is
    # checked after yt and xt have been decremented. This must be accounted
    # for.
    start_point = xt + 1, yt + 1

    # Then, we walk down, counting the steps until we encounter the lower zero
    # or edge.
    counter = 0
    yt += 1
    xt += 1
    check_bounds = lambda shape, x, y: x < shape[0] and y < shape[1] and x >= 0 and y >= 0
    while check_bounds(M.shape, xt, yt) and M[yt, xt]:
        yt += 1
        xt += 1
        counter += 1

        y_max, x_max = M.shape
        if yt >= y_max or xt >= x_max:
            break

    end_point = xt - 1, yt - 1

    return start_point, end_point, counter

from utilities import read_fasta_file

if __name__ == "__main__":
    fasta = read_fasta_file("rosalind_lcsm.txt")
    sequences = [seq for _, seq in fasta]
    s, t = sequences[0:2]

    # s = list("actgatgcaag")
    # t = list("ctgatgcata")
    # s = list("hactgu")
    # t = list("dactgo")
    # s = list("""Het besluit om voorlopig te stoppen komt nadat de redactie opnieuw een waarschuwing kreeg van de Russische mediatoezichthouder Roskomnadzor. Een recent aangenomen wet bepaalt dat alle kritiek op de oorlog in Oekraïne of de Russische overheid officieel verboden is. Er staan alleen al hoge straffen op het noemen van de woorden 'oorlog' en 'invasie'. Dat maakt onafhankelijke journalistiek in Rusland nagenoeg onmogelijk.""")
    # t = list("""Dat de krant ermee zou ophouden - in ieder geval tijdelijk - zat eraan te komen, zegt journalist Laura Starink van het kennisplatform Raam op Rusland. 'Tot voor kort was president Poetin van mening dat de oppositiepers niets voorstelde en dat hij daar niks van te vrezen had. Novaja Gazeta is een zeer kritische maar kleine krant, met een klein bereik in Rusland. Maar na de opstand in Belarus (in 2020 en 2021, red.) heeft Poetin zijn standpunt herzien.'""")
    # s = 'GTCTGGTTGTCTGAGTCCCCTACGTAAGTAGTCATTATATTAGACATGAATTAGCATCAAGCATATTGATCCTTACTTCTTTTCTAACAGCCCGGTGCCGTGCTTGTAGAAATTGATCCTCCTACCTCTCACGAAGTGCCTAGTGGTCAGGGTTATTACAAAGAATGAGAGCCGATCTTATTTCCTTGCCTGTCCTGTCATAGAGACTCGATCAACTGTGCGCGTTTGCAACTCAACGACCGGTTCCTTGCGGCGACTGAGA'

    M = np.zeros((len(s), len(t)), dtype='int')

    for i, a in enumerate(s):
        for j, b in enumerate(t):
            if a == b: M[i, j] = 1

    print(M)

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax1 = fig.add_subplot(221)

    ax1.imshow(M, interpolation='nearest')

    print("Pass", 1)
    M = addition_pass(M)
    N = M.copy()

    ax2 = fig.add_subplot(222)
    ax2.imshow(M, interpolation='nearest')


    print("Pass", 2)
    M = addition_pass(M)

    ax3 = fig.add_subplot(223)
    ax3.imshow(M, interpolation='nearest')

    # After 2 passes, the number of ones in the matrix should be stable. 
    # (Check this assumption, though!)
    # We can remove any ones if there are larger numbers in the array than one,
    # because that would indicate larger diagonals than those with a length of
    # two exist. Two-long diagonals are characterized by the ones. By removing
    # the ones, we remove the two-long diagonals.
    M[M == 1] = 0

    print("Ones removed")
    print(M)

    ax4 = fig.add_subplot(224)
    ax4.imshow(M, interpolation='nearest')

    #plt.show()

    # Continue doing passes until M.max() has only one occurrence, or until the
    # number of M.max() occurrences becomes stable in the event of two
    # diagonals of equal length.
    n_max = 2 # TODO: UUUUGGGGGHHHHHHH I hate ITM code like this :(
    prev_n_max = 0
    passes = 3
    shortenings = 0
    while n_max > 1 and prev_n_max != n_max:
        print("Pass", passes)
        M = addition_pass(M)
        if shortenings < 3:
            M = shortening_pass(M)
            shortenings += 1
        M = normalize_pass(M)

        print(M)
        print("  max:", M.max())
        print("n_max:", np.count_nonzero(M == M.max()))
        print(" ones:", np.count_nonzero(M == 1))
        print("shape:", M.shape)


        prev_n_max = n_max
        n_max = np.count_nonzero(M == M.max())
        passes += 1

    ax5 = fig.add_subplot(221)
    ax5.imshow(M, interpolation='nearest')

    plt.show()

    # Now, the highest value(s) in M indicate where the longest substrings are
    # located. Below is a plan of action for the coming sections:
    # (1) First, we get the position(s) of the maximum value(s).
    # (2) Then, we traverse up and down from these positions to determine the
    #     start point, end point, and length.
    # (3) These positions can be used to actually get the subsequence from one 
    #     of the starting sequences. 
    # (4) With this sequence in hand, we can check it against all sequences in 
    #     our haystack to see whether it occurs there.
    #     If it doesn't, we use the sequence positions to set the values on 
    #     that diagonal to 0. This allows us to apply the same routine to find 
    #     M.max() and its diagonal. Repeat from step (1) again.
    #
    # NOTE: we probably only need to do this for one, but I am implementing it
    # to do it for all of the max points, in order to verivy my assumption that
    # the same max value only occurs for the same length of substring.

    # (1) Get positions of maximum values.
    ys, xs = np.asarray(M == M.max()).nonzero()
    results = []
    # for i in range(0, n_max):
    #     y, x = ys[i], xs[i]
    for y, x in zip(ys, xs):
        # (2) For each maximum value position, walk up and down its diagonal
        #     and determine its start point, end point, and length.
        #     However, if the point lies on an already known diagonal, don't
        #     bother walking it.
        lf = lambda v: (x in range(v[0][0], v[1][0] + 1)) and (x in range([1][0], v[1][1] + 1))
        if len([r for r in results if lf(r)]) == 0:
            res = walk_up_down(N, (x, y))
            start_point, end_point, diagonal_len = res
            results.append(res)


            if diagonal_len <= 4: continue

            print(f"diagonal through {x, y}: {start_point}, {end_point}, {diagonal_len}")
            print("this corresponds to the following sequences (also they should be exactly equal!!)")
            print("s:", "".join(s[start_point[1]:end_point[1] + 1]))
            print("t:", "".join(t[start_point[0]:end_point[0] + 1]))

    # (3) Get the subsequence using these positions.
    longest_substrings = ["".join(s[start_point[1]:end_point[1]]) 
                          for start_point, end_point, _ 
                          in results]

    # (4) Go through all sequences in the file to see whether this sequence is
    #     present. Fingers crossed...
    longest_substrings.sort()
    print(longest_substrings[-1])
    max_length = 0
    hits = set()
    for ls in longest_substrings[::-1]:
        if len(ls) > max_length:
            max_length = len(ls)
            for seq in sequences:
                res = seq.find(ls)
                if res < 0:
                    break
        else:
            hits.add(ls)

    print(list(hits))
    print(max([(n, len(n)) for n in hits], key=lambda e: e[1]))

    # for id, seq in fasta:
        # r = seq.find('GTCTGGTTGTCTGAGTCCCCTACGTAAGTAGTCATTATATTAGACATGAATTAGCATCAAGCATATTGATCCTTACTTCTTTTCTAACAGCCCGGTGCCGTGCTTGTAGAAATTGATCCTCCTACCTCTCACGAAGTGCCTAGTGGTCAGGGTTATTACAAAGAATGAGAGCCGATCTTATTTCCTTGCCTGTCCTGTCATAGAGACTCGATCAACTGTGCGCGTTTGCAACTCAACGACCGGTTCCTTGCGGCGACTGAGA')
        # print(f"{id} ({r})")
