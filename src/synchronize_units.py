from copy import copy
from collections import Counter
from dataclasses import dataclass
from logzero import logger
import numpy as np
from BITS.clustering.seq import ClusteringSeq
from BITS.seq.align import EdlibRunner
from BITS.util.io import load_pickle, save_pickle
from BITS.util.proc import NoDaemonPool
from .types import TRUnit


@dataclass(eq=False)
class SyncPhase:
    """Class for synchronizing phase of the units for each read."""
    centromere_reads_fname: str = "centromere_reads.pkl"
    out_fname: str = "centromere_reads_sync.pkl"
    n_core: int = 10

    def __post_init__(self):
        self.centromere_reads = load_pickle(self.centromere_reads_fname)

    def run(self):
        with NoDaemonPool(self.n_core) as pool:
            sync_reads = pool.map(sync_units, self.centromere_reads)

        save_pickle(sync_reads, self.out_fname)


def calc_repr_units(read, ward_threshold):
    """Calculate representative units using hierarchical clustering."""
    c = ClusteringSeq([read.seq[unit.start:unit.end] for unit in read.units],
                      revcomp=False, cyclic=True if not read.synchronized else False)
    c.calc_dist_mat()
    c.cluster_hierarchical(threshold=ward_threshold)
    c.generate_consensus()
    return {df["cluster_id"]: df["sequence"] for i, df in c.cons_seqs.iterrows()}


def calc_sync_units(read, map_threshold):
    """Compute synchronized units by mapping the representative units to the read iteratively."""
    er = EdlibRunner("glocal", revcomp=False, cyclic=False)
    sync_units = []
    read_seq = copy(read.seq)
    while True:
        mappings = [(er.align(repr_unit, read_seq), repr_id)
                    for repr_id, repr_unit in sorted(read.repr_units.items())]
        diffs = [mapping.diff for mapping, repr_id in mappings]
        if np.min(diffs) >= map_threshold:
            break
        mapping, repr_id = mappings[np.argmin(diffs)]

        flatten_cigar = mapping.cigar.flatten().string
        logger.debug(mapping, repr_id)
        logger.debug(flatten_cigar)

        # Change all 'I' (= gap of a unit) at the boundaries to 'X' so that variants can be captured
        start, end = mapping.t_start, mapping.t_end
        assert flatten_cigar[0] != 'D' and flatten_cigar[-1] != 'D', "Boundary deletion happened"
        insert_len = 0   # start side
        while flatten_cigar[insert_len] == 'I':
            insert_len += 1
        while insert_len > 0 and start > 0:
            start -= 1
            insert_len -= 1
        insert_len = 0   # end side
        while flatten_cigar[-1 - insert_len] == 'I':
            insert_len += 1
        while insert_len > 0 and end < read.length:
            end += 1
            insert_len -= 1

        # TODO: remove boundary units by checking the partiality of the unit

        sync_units.append(TRUnit(start, end, id=repr_id))

        # Mask middle half sequence of the mapped region
        left = int(start + (end - start) / 4)
        right = int(end - (end - start) / 4)
        read_seq = read_seq[:left] + ('N' * (right - left)) + read_seq[right:]

    sync_units.sort(key=lambda x: x.start)

    # Resolve the conflict on the overlapping mapped regions
    er = EdlibRunner("global", revcomp=False, cyclic=False)
    for i in range(len(sync_units) - 1):
        j = i + 1
        if sync_units[i].end > sync_units[j].start:
            overlap_len = sync_units[i].end - sync_units[j].start
            logger.debug(f"conflict {sync_units[i]} and {sync_units[j]} ({overlap_len} bp)")

            # Cut out the overlapping sequeces from both units
            unit_seq = read.seq[sync_units[i].start:sync_units[i].end]
            alignment = er.align(unit_seq, read.repr_units[sync_units[i].id])
            up_seq = alignment.cigar.flatten().string[-overlap_len:]   # from unit of upper side
            
            unit_seq = read.seq[sync_units[j].start:sync_units[j].end]
            alignment = er.align(unit_seq, read.repr_units[sync_units[j].id])
            down_seq = alignment.cigar.flatten().string[:overlap_len]   # from unit of down side

            logger.debug(up_seq)
            logger.debug(down_seq)

            # Calculate the position where the total number of matches is maximized
            max_pos = 0
            max_score = 0
            for x in range(overlap_len + 1):
                score = Counter(up_seq[:x])['='] + Counter(down_seq[x:])['=']
                logger.debug(f"pos {x}, score {score}")
                if (max_score < score) or (max_score == score and up_seq[x - 1] > down_seq[x - 1]):
                    max_score = score
                    max_pos = x
            logger.debug(f"max pos {max_pos}, max score {max_score}")

            # Redefine the boundaries as the best position
            sync_units[i].end -= (overlap_len - x)
            sync_units[j].start += x

    return sync_units


def sync_units(read, ward_threshold=0.15, map_threshold=0.1):
    """Given TRRead object, synchronize the units inside it.
    <ward_threshold> is for generating representative units. 0.75 for CLR and 0.15 for CCS are recommended.
    <map_threshold> is for mapping of the representative units. 0.3 for CLR and 0.1 for CCS are recommended.
    """
    logger.info(f"read {read.id}")
    if len(read.units) <= 1:
        logger.info(f"Only <= 1 unit")
        return read

    # Compute representative units
    read.repr_units = calc_repr_units(read, ward_threshold)

    # Find all intervals similar to one of the representative units in the read
    read.units = calc_sync_units(read, map_threshold)
    read.synchronized = True

    return read
