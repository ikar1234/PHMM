"""
Parsing module to get the frquency matrix (profile), either by constructing it from a multiple
alignment or importing it as a file.
 TODO: class or function?
"""
import numpy as np
import re
from HMM.Records import Record
from typing import List
from collections import Counter
from operator import itemgetter


# TODO: other formats?
def create_pwm(file: List[str], ftype: str, alph: str) -> Record:
    """
    Create Position Weight Matrix from Multiple Alignment
    :param file: alignment
    :param ftype: format of MSA
    :param alph: alphabet
    """
    # TODO: many alingments in one file?
    header = file[0].strip().split(' ')
    n_seq = int(header[0])
    # extract ids
    ids = [re.findall('.*? ', file[i])[0] for i in range(1, n_seq + 1)]
    # iterate through every row i modulo number of sequences and extract the sequence
    # then remove spaces
    # do this for every sequence k from 1 to n_seq
    seq = [re.sub(' ', '', ''.join([re.findall(' (.*)', file[i])[0].strip() for i in range(k, len(file), n_seq + 1)]))
           for k in range(1, n_seq + 1)]
    # TODO: currently, meta data is parsed into MEME format; change?
    # sequence length
    seq_len = len(seq[0])
    if alph.upper().startswith('A'):
        alength = 4
        al = "A C G T".split(' ')
    else:
        # add ? character
        alength = 21
        al = "A C D E F G H I K L M N P Q R S T V W Y ?".split(' ')
    matrix = np.zeros((alength, seq_len))
    for i in range(seq_len):
        c = Counter([s[i] for s in seq])
        # add letters which occur 0 times
        for a in al:
            if c.get(a, 0) == 0:
                c[a] = 0
        c = sorted(c.items(), key=itemgetter(0))
        # fill matrix column-wise
        matrix[:, i] = np.array([x[1] for x in c])/n_seq
    # add number of sequences, alphabet length and sequence length
    meta = {'alength': alength, 'nsites': n_seq, 'E': 0}
    # use first sequence for name
    name = ids[0]
    record = Record(name=name, matrix=matrix, ftype=ftype, meta=meta)
    return record


# TODO: some matrices come normalized, other unnormalized
def parse(filepath, ftype="fasta", alph='dna') -> Record:
    """
    :param filepath: Path to file
    :param ftype: File type
    :return: Name of motif/etc. and an un/normalized count matrix
    """
    if ftype.upper() not in ['PHYLIP', 'JASPAR', 'MEME']:
        print(f"Input should be one of the following types: PHYLIP, ")
        record = Record()
        return record

    with open(filepath, "r") as f:
        file = f.readlines()

    if ftype.upper() == "MEME":
        try:
            name = [line for line in file if line.startswith("MOTIF")][0]
        except ValueError:
            raise ValueError("Cannot find motif name.")
        name = ' '.join(name.split()[1:])
        start = [x for x in range(len(file)) if file[x].startswith("letter-probability matrix")][0]
        # capture meta information
        reg = re.findall(pattern=r"\w+= \d+",
                         string=file[start])
        meta = {r.split('=')[0]: int(r.split('=')[1]) for r in reg}
        matrix = []
        # if meta info does not contain sequence length, go to the end of the file
        for i in range(start + 1, start + 1 + meta.get('w', len(file) - start - 1)):
            if file[i].startswith(r"\w+"):
                break
            nums = file[i].strip().split()
            matrix.append([float(x) for x in nums])
        # take the transpose, since the sequence should be 'horizontal'
        matrix = np.array(matrix).transpose()
        record = Record(name, matrix, ftype="MEME", meta=meta)
        return record

    elif ftype.upper() == "JASPAR":
        # TODO: maybe more?
        name = file[0][1:].strip()
        matrix = []
        for i in range(1, len(file)):
            # print(file[i])
            reg = re.search(pattern=r"(\d+\s+)+",
                            string=file[i])
            nums = [float(x) for x in reg.group(0).split()]
            matrix.append(nums)
        matrix = np.array(matrix)

        record = Record(name, matrix, ftype="JASPAR")
    # else: file is in PHYLIP format
    else:
        record = create_pwm(file, ftype=ftype.upper(), alph=alph)

    return record


p = "C:\\Users\\user\\PycharmProjects\\Helloadvanced\\HMM\\phylip.txt"

parse(p, ftype="phylip")
