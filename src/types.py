from dataclasses import dataclass, astuple
from typing import List, Dict
import numpy as np
from BITS.seq.utils import revcomp


@dataclass(eq=False)
class Read:
    """Class for a read."""
    seq : str
    id  : int = None   # for DAZZ_DB (read) or cluster ID (TR representative unit)
    name: str = None   # for fasta

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

    instance variables:
      @ repr_id <repr_id> [None] : ID of representative unit which this unit belongs to
      @ strand  <int>     [None] : 0 (forward) or 1 (revcomp)

    <strand>(<TRRead.seq>[start:end]) aligns to <TRRead.repr_units>[<repr_id>].
    Representative units and TR intervals are always defined in the forward direction.
    """
    repr_id: int = None
    strand : int = None


@dataclass(eq=False)
class TRRead(Read):
    """Class for a read with TRs. Multiple TRs in a read are not distinguished here.

    instance variables:
      @ alignments   <List[SelfAlignment]> [None]  : outout of datander
      @ trs          <List[ReadInterval]>  [None]  : output of datander
      @ units        <List[TRUnit]>        [None]  : initially computed by datruf
      @ repr_units   <Dict[int, str]>      [None]  : {repr_id: sequence}; only forward sequence is OK
      @ synchronized <bool>                [False] : whether or not <units> are
      @ quals        <np.ndarray>          [None]  : Positional QVs
    """
    alignments  : List[SelfAlignment] = None
    trs         : List[ReadInterval]  = None
    units       : List[TRUnit]        = None
    repr_units  : Dict[int, str]      = None
    synchronized: bool                = False
    quals       : np.ndarray          = None

    def _check_unit_strands(self):
        """Check if all the unit strands are defined."""
        for unit in self.units:
            assert unit.strand is not None, "Unit strand is not defined."

    @property
    def unit_seqs(self, forward=False):
        """Return TR unit sequences. If <forward> is True, the orientations of the units are modified 
        so that these are same as those of <repr_units>."""
        if forward:
            self._check_unit_strands()

        return [self.seq[unit.start:unit.end] if not forward or unit.strand == 0
                else revcomp(self.seq[unit.start:unit.end])
                for unit in self.units]

    @property
    def unit_quals(self, forward=False):
        if forward:
            self._check_unit_strands()

        return [self.quals[unit.start:unit.end] if not forward or unit.strand == 0
                else np.flip(self.quals[unit.start:unit.end])
                for unit in self.units]


@dataclass
class Overlap:
    """Class for an overlap between two reads.

    positional arguments:
      @ a_read_id <int>   : `a_read[a_start:a_end]` is the overlapping sequence.
      @ b_read_id <int>   : `strand(b_read)[b_start:b_end]` is the overlapping sequence.
      @ strand    <int>   : Must be 0 or 1
      @ a_start   <int>   : Start position of the overlap on `a_read`. 0-indexed.
      @ a_end     <int>   : End position on `a_read`.
      @ a_len     <int>   : Length of the `a_read`.
      @ b_start   <int>   : NOTE: `[a|b]_[start|end]` are always defined for FORWARD overlapping sequences.
      @ b_end     <int>
      @ b_len     <int>
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

    def astuple(self):
        return astuple(self)
