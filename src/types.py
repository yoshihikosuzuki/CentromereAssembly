from dataclasses import dataclass, field
from typing import List, Dict
import numpy as np


@dataclass(frozen=True)
class SelfAlignment:
    abpos: int
    aepos: int
    bbpos: int
    bepos: int
        
    @property
    def distance(self):
        return self.abpos - self.bbpos
        
    @property
    def slope(self):
        return round((self.aepos - self.abpos) / (self.bepos - self.bbpos), 3)


@dataclass(eq=False)
class TRUnit:
    """Normally used as an instance in a Read or TR object."""
    start : int
    end   : int
    type  : str = None   # "complete" or "partial"
    id    : int = None   # for clustering of units   # TODO: change name based on the clustering method

    @property
    def length(self):
        return self.end - self.start


@dataclass(eq=False)
class TR:
    """Normally used as an instance in a Read object."""
    start: int
    end: int
    units: List[TRUnit] = None

    @property
    def cv_ulen(self):
        ulens = [unit.length for unit in self.units]
        return round(np.std(ulens, ddof=1) / np.mean(ulens), 3)


@dataclass(eq=False)
class Read:
    seq  : str = None   # can be None (abstract read)
    id   : int = None   # for DAZZ_DB
    name : str = None   # for fasta


@dataclass(eq=False)
class TRRead(Read):
    alignments : List[SelfAlignment] = None   # used in datruf
    trs        : List[TR]            = None
    repr_units : Dict[int, str]      = None   # {cluster_id: str}
    sync_units : List[TRUnit]        = None   # synchronized/classified units
