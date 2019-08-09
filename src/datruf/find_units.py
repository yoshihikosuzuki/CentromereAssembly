from multiprocessing import Pool
from interval import interval
from logzero import logger
from BITS.util.io import save_pickle
from BITS.util.interval import intvl_len, subtract_intvl
from .io import load_dumps, load_paths
from ..types import TR, TRUnit, TRRead


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
    """Call <find_units_single> for each read whose id is in [<start_dbid>:<end_dbid> + 1]."""
    # Load TR reads with data of TR intervals and all self alignments
    tr_reads_dump = load_dumps(start_dbid, end_dbid, db_fname, las_fname)

    # For each read, calculate the unit intervals
    tr_reads = []
    for tr_read_dump in tr_reads_dump:
        tr_read = find_units_single(tr_read_dump, db_fname, las_fname)
        if tr_read is not None:
            tr_reads.append(tr_read)
    return tr_reads


def find_units_single(tr_read_dump, db_fname, las_fname, max_cv=0.1):
    """Core function of datruf.
    Find the best set of self alignments and split the TR intervals induced by the alignments into units."""
    trs = []
    inner_alignments = find_inner_alignments(tr_read_dump)
    inner_paths = load_paths(tr_read_dump, inner_alignments, db_fname, las_fname)
    for alignment, fcigar in sorted(inner_paths.items(), key=lambda x: x[0].bbpos):
        tr_units = split_tr(alignment.abpos, alignment.bbpos, fcigar)
        if len(tr_units) == 1:   # at least duplication is required
            continue

        tr = TR(start=alignment.bbpos, end=alignment.aepos, units=tr_units)
        if tr.cv_ulen >= max_cv:   # probably unit length is too short
            continue
        
        trs.append(tr)

    return TRRead(id=tr_read_dump.id, trs=trs) if len(trs) > 0 else None


def find_inner_alignments(tr_read_dump, min_len=1000):
    """Extract a set of non-overlapping most inner self alignments.
    <min_len> defines the required overlap length with yet uncovered TR region."""
    uncovered = interval(*tr_read_dump.trs)
    inner_alignments = set()
    for alignment in tr_read_dump.alignments:   # in order of distance
        if intvl_len(uncovered) < min_len:
            break
        intersect = uncovered & interval[alignment.bbpos, alignment.aepos]
        uncovered = subtract_intvl(uncovered, interval[alignment.bbpos, alignment.aepos])
        if (intvl_len(intersect) >= min_len
            and 0.95 <= alignment.slope <= 1.05   # eliminate abnornal slope
            and alignment.abpos <= alignment.bepos):   # at least duplication
            inner_alignments.add(alignment)   # TODO: add only intersection is better?
    logger.debug(f"inners: {inner_alignments}")
    return inner_alignments


def split_tr(abpos, bbpos, fcigar):
    # Split TR interval into unit intervals given <fcigar> specifying self alignment
    # <abpos> corresponds to the start position of the first unit
    # <bbpos> does the second
    apos, bpos = abpos, bbpos
    tr_units = [TRUnit(start=bpos, end=apos)]
    # Iteratively find max{ax} such that bx == (last unit end)
    for i, c in enumerate(fcigar):
        if c != 'I':
            apos += 1
        if c != 'D':
            bpos += 1
        if bpos == tr_units[-1][1] and (i == len(fcigar) - 1 or fcigar[i + 1] != 'D'):
            tr_units.append(TRUnit(start=tr_units[-1][1], end=apos))
    return tr_units
