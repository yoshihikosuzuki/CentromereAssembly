from multiprocessing import Pool
import numpy as np
from interval import interval
from logzero import logger
from BITS.util.io import save_pickle
from BITS.util.interval import intvl_len, subtract_intvl
from .io import load_tr_reads, load_paths
from ..types import TRUnit


def find_units(start_dbid, end_dbid, n_core, db_fname, las_fname, out_fname):
    """Split the tasks ranging from <start_dbid> and <end_dbid> into <n_core> sub-tasks.
    Each splitted sub-task is then parallely executed, and output into a file.
    A function named "find_units_single" above is actually the core function.
    """
    if n_core == 1:
        tr_reads = find_units_multi(start_dbid, end_dbid, db_fname, las_fname)
    else:
        unit_n = -(-(end_dbid - start_dbid + 1) // n_core)
        args = [(start_dbid + i * unit_n,   # start_dbid
                 min([start_dbid + (i + 1) * unit_n - 1, end_dbid]),   # end_dbid
                 db_fname, las_fname)
                for i in range(n_core)]
        tr_reads = []
        with Pool(n_core) as pool:
            for tr_reads_unit in pool.starmap(find_units_multi, args):
                tr_reads += tr_reads_unit

    save_pickle(tr_reads, out_fname)


def find_units_multi(start_dbid, end_dbid, db_fname, las_fname):
    """Call <find_units_single> for each read whose id is in [<start_dbid>:<end_dbid> + 1].
    This returns all the TR reads even when CV of the unit lengths is large although units are not determined.
    """
    # Load TR reads with data of TR intervals and all self alignments
    reads = load_tr_reads(start_dbid, end_dbid, db_fname, las_fname)

    # For each read, calculate the unit intervals
    # NOTE: <read.units> can be empty list (i.e. TRs are detected but alignments are noisy or too short)
    for read in reads:
        read.units = find_units_single(read, db_fname, las_fname)

    return reads


def find_units_single(read, db_fname, las_fname, max_cv=0.1):
    """Core function of datruf.
    Find the best set of self alignments and split the TR intervals induced by the alignments into units.
    """
    all_units = []
    # Determine a set of self alignments from which units are cut out
    inner_alignments = find_inner_alignments(read)
    # Load flattten CIGAR strings of the selected alignments
    inner_paths = load_paths(read, inner_alignments, db_fname, las_fname)

    for alignment, fcigar in inner_paths.items():
        # Compute unit intervals based on the reflecting snake
        # between the read and the self alignment
        units = split_tr(alignment.ab, alignment.bb, fcigar)
        if len(units) == 1:   # at least duplication is required
            continue

        # Exclude TRs with abnormal CV (probably due to short unit length)
        # and then add the units
        ulens = [unit.length for unit in units]
        cv_ulen = round(np.std(ulens, ddof=1) / np.mean(ulens), 3)
        if cv_ulen >= max_cv:
            continue
        all_units += units

    # TODO: remove "contained units"
    return all_units


def find_inner_alignments(read, min_len=1000):
    """Extract a set of non-overlapping most inner self alignments.
    <min_len> defines the required overlap length with yet uncovered TR region."""
    uncovered = interval(*[(tr.start, tr.end) for tr in read.trs])
    inner_alignments = set()
    for alignment in read.alignments:   # in order of distance
        if intvl_len(uncovered) < min_len:
            break
        intersect = uncovered & interval[alignment.bb, alignment.ae]
        uncovered = subtract_intvl(uncovered, interval[alignment.bb, alignment.ae])
        if (intvl_len(intersect) >= min_len
            and 0.95 <= alignment.slope <= 1.05   # eliminate abnornal slope
            and alignment.ab <= alignment.be):   # at least duplication
            inner_alignments.add(alignment)   # TODO: add only intersection is better?
    logger.debug(f"inners: {inner_alignments}")
    return inner_alignments


def split_tr(ab, bb, fcigar):
    """Split TR interval into unit intervals given <fcigar> specifying self alignment
    <ab> corresponds to the start position of the first unit
    <bb> does the second"""
    apos, bpos = ab, bb
    tr_units = [TRUnit(start=bpos, end=apos)]
    # Iteratively find max{ax} such that bx == (last unit end)
    for i, c in enumerate(fcigar):
        if c != 'I':
            apos += 1
        if c != 'D':
            bpos += 1
        if bpos == tr_units[-1].end and (i == len(fcigar) - 1 or fcigar[i + 1] != 'D'):
            tr_units.append(TRUnit(start=tr_units[-1].end, end=apos))
    return tr_units
