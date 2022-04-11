from urllib import request
from urllib.error import HTTPError, URLError
from collections import namedtuple

FastaEntry = namedtuple('FastaEntry', 'id seq')

def _str_str_to_FastaEntry(entry: list[str]) -> FastaEntry:
    id, seq, *_ = entry
    return FastaEntry(id, seq)

def get_uniprot_fasta(uniprot_id: str) -> FastaEntry:
    # TODO: use a proper url constructer for this.
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"

    with request.urlopen(url) as response:
        body = response.read()

    # We ignore any other entries, because we can assume only one fasta entry 
    # is present in the response from UniProt.
    entry, *_ = read_fasta(body.decode('utf-8'))
    return entry

def get_rcsb_pdb_fasta(pdb_id: str) -> list[FastaEntry]:
    # TODO: use a proper url constructer for this.
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

    with request.urlopen(url) as response:
        body = response.read()

    # Responses from RCSB PDB can contain a number of entries, each
    # representing the sequence of a different chain. We return each of the
    # entries, rather than a single FastaEntry as get_uniprot_fasta() does.
    entries = read_fasta(body.decode('utf-8'))
    return entries

def file_to_string(path: str) -> str:
    data = ""
    with open(path, "r") as f:
        data = f.read()

    return data

def string_to_file(s: str, path: str):
    with open(path, "w") as f:
        f.write(s)

def read_fasta(fasta: str) -> list[FastaEntry]:
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

    entries = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            # The line is the start of a new fasta entry.
            ident = line[1:].strip()
            entries.append([ident, ""])
        else:
            entries[-1][1] += line.strip()

    return [FastaEntry(id, seq) for id, seq in entries]

def read_fasta_file(path: str) -> list[FastaEntry]:
    fasta = file_to_string(path)
    return read_fasta(fasta)
