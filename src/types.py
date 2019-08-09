from dataclasses import dataclass, field
from typing import List, Dict
import numpy as np


@dataclass(frozen=True)
class SelfAlignment:
    """Class for a self alignment. Used in datruf."""
    ab: int
    ae: int
    bb: int
    be: int
        
    @property
    def distance(self):
        return self.ab - self.bb
        
    @property
    def slope(self):
        return round((self.ae - self.ab) / (self.be - self.bb), 3)


@dataclass(eq=False)
class ReadInterval:
    """Class for an abstract interval of a read.
    A tandem repeat (not unit) is represented using this."""
    start: int
    end: int

    @property
    def length(self):
        return self.end - self.start


@dataclass(eq=False)
class TRUnit(ReadInterval):
    """Class for a tandem repeat unit. Equal to ReadInterval with some properties.
    Normally used as instance variable of TRRead."""
    complete : bool = True
    id       : int  = None   # for clustering of units   # TODO: change name based on the clustering method


@dataclass(eq=False)
class Read:
    """Class for a read."""
    seq  : str = None   # can be None (abstract read)
    id   : int = None   # for DAZZ_DB
    name : str = None   # for fasta


@dataclass(repr=False, eq=False)
class ReadDump(Read):
    """Class for a read plus dump data. Used in datruf."""
    trs: List[ReadInterval]
    alignments: List[SelfAlignment]


@dataclass(eq=False)
class TRRead(Read):
    """Class for a read with TRs. Multiple TRs in a read are not distinguished here."""
    units      : List[TRUnit]   = None
    repr_units : Dict[int, str] = None   # {cluster_id: str}
