"""
Align sequences from two fasta files and present plots.

Usage:
    aligner [options] [alignment] [fasta file] [fasta file]

Available options are:
    -h, --help          Show this help
    --plot-alignment    Plot the alignment matrix
    --plot-filtering    Plot the filtering matrix
    --substrings        Report 10 longest substrings

Available alignment options:
    hamming             Apply Hamming distance alignment
    windowed-blosum     Apply a windowed BLOSUM62 alignment
    blosum              Apply a simple pairwise BLOSUM62 alignment
    single              Apply a simple pairwise identity alignment
    hydropathy          Apply an alignment based on Kyte-Doolittle hydropathy

Contact:
    https://github.com/koenwestendorp/aligner

Version:
    aligner v0.1.5
"""

import sys

from aligner import aligner, utilities
import matplotlib.pyplot as plt

def main():
    """
    Align sequences from two fasta files and present plots.
    """

    # Rather elegant method for quickly parsing simple command line arguments I
    # found: https://github.com/realpython/reader/blob/master/reader/__main__.py
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

    if "-h" in opts or "--help" in opts:
        print(__doc__)
        raise SystemExit()

    # Check whether the number of arguments is 2.
    if len(args) == 3:
        alignments = ["hamming", "windowed-blosum", "blosum", "single", 
                      "hydropathy"]

        alignment, *fasta_files = args
        if alignment not in alignments:
            print(alignment, "is not a known alignment.")
            print("See --help for a list of alignment methods and usage.")
            raise SystemExit()

        # Assuming single-entry fasta files, we take the first entry.
        entries = []
        for file_path in fasta_files:
            print(f"Reading '{file_path}' ...")
            entry = utilities.read_fasta_file(file_path)[0]
            entries.append(entry)

        a = aligner.Aligner(*entries, window_size=3, threshold=3)

        # Construct the alignment matrix according to the `alignment` argument.
        if alignment == "hamming":
            a.construct_windowed_alignment_matrix(scoring_func=aligner.hamming_distance_inverse)
        if alignment == "windowed-blosum":
            a.construct_windowed_alignment_matrix(scoring_func=aligner.blosum_62_windowed)
        if alignment == "blosum":
            a.construct_alignment_matrix(scoring_func=aligner.blosum_62)
        if alignment == "single":
            a.construct_alignment_matrix()
        if alignment == "hydropathy":
            a.construct_alignment_matrix(scoring_func=aligner.hydrophobicity_score)  # type: ignore

        a.addition_pass()
        a.addition_pass()

        if "--plot-alignment" in opts:
            a.plot_alignment_matrix()
        if "--plot-filtering" in opts:
            a.plot_filter_matrix()
        if "--substrings" in opts:
            a.long_substrings_report()
    else:
        print("Invalid usage. Try")
        print("\taligner [options] [alignment] [fasta file] [fasta file]")
        print("or see --help for usage.")
        raise SystemExit()

if __name__ == "__main__":
    main()
