from logzero import logger
from BITS.seq.align import EdlibRunner
from BITS.seq.utils import revcomp
from BITS.util.io import save_pickle, load_pickle
from BITS.util.proc import run_command, NoDaemonPool


def read_to_boundary_k_monomers(read, k=2, offset=1):
    """Extract <k>-monomers (including sequences between monomers) of the boundaries.
    Sequences of both orientations are returned. Therefore, the total number of <k>-monomers are 4.
    <offset> boundary monomers are skipped to avoid noise on sequence."""
    if len(read.units) < k + 2 * offset:
        logger.info(f"Read {read.id}: # of units not sufficient. Skip.")
        return []
    
    k_units = [read.seq[read.units[offset].start:read.units[offset + k - 1].end],   # prefix boundary
               read.seq[read.units[-offset - k].start:read.units[-offset - 1].end]]   # suffix boundary
    k_units += [revcomp(k_unit) for k_unit in k_units]   # add revcomp seqs

    return k_units


def reads_to_all_kplus1_monomers(reads, k=2):
    """List up all (<k>+1)-monomers (including sequences between monomers) contained in <reads>.
    Only forward sequences are returned."""
    all_kplus1_units = []
    for read in reads:
        if len(read.units) < k + 1:
            logger.info(f"Read {read.id}: # of units not sufficient. Skip.")
            continue
        all_kplus1_units += [(read.seq[read.units[i].start:read.units[i + k].end], read.id)
                             for i in range(len(read.units) - k)]
    return all_kplus1_units


def sva_k_monomer_ovlp(read, all_kplus1_units, er, max_diff=0.015):
    """Given a single read, map 4 boundary k-units of the read to all the (k + 1)-units of all the reads
    and return read pairs with good mapping."""
    ovlps = set()
    query_read_id = read.id
    for boundary_k_unit in read_to_boundary_k_monomers(read):
        for kplus1_unit, target_read_id in all_kplus1_units:
            if query_read_id == target_read_id:
                continue
            if tuple(sorted([query_read_id, target_read_id])) in ovlps:
                continue
            alignment = er.align(boundary_k_unit, kplus1_unit)
            if alignment.diff < max_diff:
                logger.info(f"{query_read_id}~{target_read_id} ({round(100 * alignment.diff, 1)}% diff)")
                ovlps.add(tuple(sorted([query_read_id, target_read_id])))
    return ovlps


if __name__ == "__main__":
    import sys
    i, n_distribute, n_core, out_fname = sys.argv[1:]
    i, n_distribute, n_core = list(map(int, [i, n_distribute, n_core]))

    centromere_reads = load_pickle("centromere_reads.pkl")
    all_kplus1_units = reads_to_all_kplus1_monomers(centromere_reads)
    er = EdlibRunner("glocal", revcomp=False, cyclic=False)

    unit_n = -(-len(centromere_reads) // n_distribute)
    centromere_reads = centromere_reads[i * unit_n:(i + 1) * unit_n]

    overlaps = set()
    with NoDaemonPool(n_core) as pool:
        for ret in pool.starmap(sva_k_monomer_ovlp, [(centromere_read, all_kplus1_units, er)
                                                     for centromere_read in centromere_reads]):
            overlaps.update(ret)

    save_pickle(overlaps, out_fname)
