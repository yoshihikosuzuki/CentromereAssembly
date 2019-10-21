from dataclasses import dataclass, astuple
from typing import List, Dict
import numpy as np
from BITS.seq.utils import revcomp_seq


@dataclass(eq=False)
class Read:
    """Class for a read.

    positional instance variables:
      @ seq <str>

    optional instance variables:
      @ id     <int> [None] : DAZZ_DB id
      @ name   <str> [None] : Fasta header
      @ strand <int> [0]    : 0 (forward) or 1 (revcomp)
    """
    seq   : str
    id    : int = None
    name  : str = None
    strand: int = 0

    def __post_init__(self):
        assert self.strand == 0 or self.strand == 1, "Strand must be 0 or 1"

    @property
    def length(self):
        return len(self.seq)


@dataclass(frozen=True)
class SelfAlignment:
    """Class for a self alignment calculated by datander. Used in datruf."""
    ab: int
    ae: int
    bb: int
    be: int

    @property
    def astuple(self):
        return astuple(self)

    @property
    def distance(self):
        """From diagonal. Same as the length of the first unit."""
        return self.ab - self.bb

    @property
    def slope(self):
        """One metric on how the self alignment is noisy. Units are unreliable if the value is large,
        although the TR interval [bb, ae] is still somewhat reliable."""
        return round((self.ae - self.ab) / (self.be - self.bb), 3)


@dataclass(eq=False)
class ReadInterval:
    """Class for an interval within a read.
    A tandem repeat (not unit) is represented by this class."""
    start: int
    end  : int

    @property
    def length(self):
        return self.end - self.start


@dataclass(eq=False)
class TRUnit(ReadInterval):
    """Class for a tandem repeat unit. This class does not store the sequence and is used as an instance
    variable of TRRead.

    positional instance variables:
      @ start <int>
      @ end   <int>

    optional instance variables:
      @ repr_id <repr_id> [None] : ID of representative unit which this unit belongs to
      @ strand  <int>     [None] : 0 (forward) or 1 (revcomp)

    `strand(read.seq[start:end])` aligns to `read.repr_units[repr_id]`.
    Representative units and TR intervals are always defined in the forward direction.
    """
    repr_id: int = None
    strand : int = None


@dataclass(eq=False)
class TRRead(Read):
    """Class for a read with TRs.

    positional instance variables:
      @ seq <str>

    optional instance variables:
      @ id           <int>                 [None]  : DAZZ_DB id
      @ name         <str>                 [None]  : Fasta header
      @ strand       <int>                 [None]  : 0 (forward) or 1 (revcomp)
      @ alignments   <List[SelfAlignment]> [None]  : Outout of datander
      @ trs          <List[ReadInterval]>  [None]  : Output of datander
      @ units        <List[TRUnit]>        [None]  : Initially computed by datruf
      @ repr_units   <Dict[int, str]>      [None]  : `{repr_id: sequence}`. Assumed only forward sequences.
      @ synchronized <bool>                [False] : Whether or not `self.units` are
      @ quals        <np.ndarray>          [None]  : Positional QVs
    """
    alignments  : List[SelfAlignment] = None
    trs         : List[ReadInterval]  = None
    units       : List[TRUnit]        = None
    repr_units  : Dict[int, str]      = None
    synchronized: bool                = False
    quals       : np.ndarray          = None

    @property
    def unit_seqs(self, forward=False):
        """Return TR unit sequences. If <forward> is True, the orientations of the units are modified 
        so that these are same as those of <repr_units>."""
        return [self.seq[unit.start:unit.end] if not forward or unit.strand == 0
                else revcomp_seq(self.seq[unit.start:unit.end])
                for unit in self.units]

    @property
    def unit_quals(self, forward=False):
        return [self.quals[unit.start:unit.end] if not forward or unit.strand == 0
                else np.flip(self.quals[unit.start:unit.end])
                for unit in self.units]


def revcomp_read(read):
    """Return reverse complement of <read> as a new object."""
    # TODO: revcomp `alignments` and `trs`
    return TRRead(seq=revcomp_seq(read.seq), id=read.id, name=read.name, strand=1 - read.strand,
                  units=[TRUnit(start=read.length - unit.end,
                                end=read.length - unit.start,
                                repr_id=unit.repr_id,
                                strand=(None if unit.strand is None else 1 - unit.strand))
                         for unit in reversed(read.units)],
                  synchronized=read.synchronized,
                  repr_units=read.repr_units,
                  quals=None if read.quals is None else np.flip(read.quals))


@dataclass
class Overlap:
    """Class for an overlap between two reads.

    positional instance variables:
      @ a_read_id <int>   : `a_read[a_start:a_end]` is the overlapping sequence.
      @ b_read_id <int>   : `strand(b_read)[b_start:b_end]` is the overlapping sequence.
      @ strand    <int>   : 0 (forward) or 1 (revcomp)
      @ a_start   <int>   : Start position of the overlap on `a_read`. 0-indexed.
      @ a_end     <int>   : End position on `a_read`.
      @ a_len     <int>   : Length of the `a_read`.
      @ b_start   <int>   : NOTE: `[a|b]_[start|end]` are always defined for FORWARD overlapping sequences.
      @ b_end     <int>   : 
      @ b_len     <int>   : 
      @ diff      <float> : Percent dissimilarity at the overlap.
    """
    a_read_id: int
    b_read_id: int
    strand   : int
    a_start  : int
    a_end    : int
    a_len    : int
    b_start  : int
    b_end    : int
    b_len    : int
    diff     : float

    def __post_init__(self):
        assert self.strand == 0 or self.strand == 1, "`strand` must be 0 or 1"

    def astuple(self):
        return astuple(self)
