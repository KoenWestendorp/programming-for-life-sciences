def file_to_string(path: str) -> str:
    data = ""
    with open(path, "r") as f:
        data = f.read()

    return data

def string_to_file(s: str, path: str):
    with open(path, "w") as f:
        f.write(s)

def read_fasta(fasta: str) -> [[str, str]]:
    '''
    Returns a list of lists countaining the identifier string and the sequence
    content.

    If the fasta does not contain any entries, an empty list is returned.

    Parameters
    ----------
    fasta : str
        A string of fasta data, which may contain multiple entries.

    Returns
    -------
    [[str, str]]
        List of [identifier, sequence].
    '''

    ret = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            # The line is the start of a new fasta entry.
            ident = line[1:].strip()
            ret.append([ident, ""])
        else:
            ret[-1][1] += line.strip()

    return ret

def read_fasta_test():
    fasta = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
    data = [
        ['Rosalind_6404', 'CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG'], 
        ['Rosalind_5959', 'CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC'], 
        ['Rosalind_0808', 'CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT'],
    ]
    assert read_fasta(fasta) == data

def read_fasta_file(path: str) -> [[str, str]]:
    fasta = file_to_string(path)
    return read_fasta(fasta)

if __name__ == "__main__":
    read_fasta_test()
