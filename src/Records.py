"""
Similar to the Record class from Bio.motifs.records, but I'll tyr to make it more flexible for my purposes
"""
import numpy as np
from typing import Dict


class Record:
    """
    Record class which holds the information about the parsed PWM.
    """
    name: str
    matrix: np.ndarray
    ftype: str
    meta: Dict
    __slots__ = ("name", "matrix", "ftype", "meta")

    def __init__(self, name=None, matrix=None, ftype=None, meta=None):
        self.name = name
        self.matrix = matrix
        self.ftype = ftype
        self.meta = meta

    def __str__(self):
        s = f"Record \n" \
            f"Name: {self.name} \n" \
            f"Position Weight Matrix:\n {self.matrix} \n" \
            f"Input format: {self.ftype} \n" \
            f"Meta data: {', '.join([k+'= '+str(v) for k,v in self.meta.items()])}"
        return s

    def todict(self) -> Dict:
        return vars(self)
