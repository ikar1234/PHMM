"""
Used for storing the sequences or the multiple alignment.
"""
from typing import Dict, List


class Record:
    """
    Record class which holds the information about the parsed PWM.
    """
    name: str
    seqs: List[str]
    alph: str
    aligned: bool
    __slots__ = ("name", 'seqs', 'alph', 'aligned')

    def __init__(self, name=None, seqs=None, alph="dna", aligned=False):
        self.name = name
        self.seqs = seqs
        self.alph = alph
        self.aligned = aligned

    def __str__(self):
        s = f"Record \n" \
            f"Name: {self.name} \n" \
            f"Alphabet: {self.alph} \n"
        return s

    def todict(self) -> Dict:
        return vars(self)
