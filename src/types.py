from dataclasses import dataclass
from typing import List


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
class SeqInterval:
    """Class for an abstract interval of a seqeunce including read.
    A tandem repeat (not unit) is represented by this class."""
    start: int
    end: int

    @property
    def length(self):
        return self.end - self.start


@dataclass(eq=False)
class TRUnit(SeqInterval):
    """Class for a tandem repeat unit. Equal to SeqInterval with some properties.
    Normally used as instance variable of TRRead."""
    repr_id : int = None   # ID of representative unit to which this TRUnit belongs


@dataclass(eq=False)
class Sequence:
    """Class for a sequence including read."""
    seq    : str
    id     : int = None   # for DAZZ_DB (read) or cluster ID (TR representative unit)
    name   : str = None   # for fasta

    @property
    def length(self):
        return len(self.seq)


@dataclass(eq=False)
class TRRead(Sequence):
    """Class for a read with TRs. Multiple TRs in a read are not distinguished here."""
    alignments   : List[SelfAlignment] = None
    trs          : List[SeqInterval]  = None
    units        : List[TRUnit]        = None
    repr_units   : List[Sequence]      = None   # only forward sequences should be registered;
                                                # that is, for the master units of all the reads you
                                                # should register both forward and reverse complement
                                                # of the master unit sequence with IDs like 0 and 1
    synchronized : bool                = False   # whether or not <units> are
