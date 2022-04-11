import aligner.aligner as al
import numpy as np

def test_blosum_62():
    assert al.blosum_62('G', 'Q') == -2, \
        "Score for GQ must be -2."

def test_blosum_62_windowed():
    assert al.blosum_62_windowed(['E', 'F', 'G'], ['H', 'I', 'K']) == -2, \
        "Score for EFG against HIK must be -2."

def test_hamming_distance():
    assert al.hamming_distance(['A', 'M', 'Y'], ['I', 'M', 'A']) == 2, \
        "Hamming distance between AMY and IMA must be 2."
    assert al.hamming_distance(['S'], ['I', 'S', 'L']) == 1, \
        "Hamming distance between S and ISL must be 1."

def test_hamming_distance_inverse():
    assert al.hamming_distance_inverse(['A', 'M', 'Y'], ['I', 'M', 'A']) == 1, \
        "Inverse Hamming distance between AMY and IMA must be 2."
    assert al.hamming_distance_inverse(['S'], ['I', 'S', 'L']) == 0, \
        "Inverse Hamming distance between S and ISL must be 0."

def test_kyte_doolittle_hydropathy_factor():
    assert al.kyte_doolittle_hydropathy_factor('G') == 4.1

def test_hydrophobicity_score():
    assert al.hydrophobicity_score('g', 'A') == 22

def test_windowed():
    assert al.windowed(['a', 'b', 'c', 'd', 'e'], 3) == [['a', 'b', 'c'], 
                                                         ['b', 'c', 'd'], 
                                                         ['c', 'd', 'e'], 
                                                         ['d', 'e'], 
                                                         ['e']]

def test_walk_up_down():
    M = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 2, 0, 0, 0, 0, 0, 4, 0, 0],
                  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    pos1 = (4, 4)
    assert al.walk_up_down(M, pos1) == ((2, 2), (5, 5), 4)
    pos2 = (9, 3)
    assert al.walk_up_down(M, pos2) == ((8, 2), (10, 4), 3)
