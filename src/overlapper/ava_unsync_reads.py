import argparse
from dataclasses import dataclass
from multiprocessing import Pool
import numpy as np
from logzero import logger
from BITS.seq.kmer import seq_to_forward_kmer_spectrum
from BITS.util.io import save_pickle, load_pickle
from BITS.util.proc import run_command
from BITS.util.scheduler import Scheduler
from ..types import revcomp_read
from .svs_unsync_reads import svs_overlap

out_dir        = "ava_unsync"
out_prefix     = "ovlps"
scatter_prefix = "run_ava_unsync"
gather_fname   = f"{out_dir}/gather.sh"
log_prefix      = f"{out_dir}/log"


@dataclass(eq=False)
class UnsyncReadsOverlapper:
    """Class for executing all-vs-all overlap between unsynchronized TR reads.

    positional arguments:
      @ n_distribute <int> : Number of jobs
      @ n_core       <int> : Number of cores per job

    optional arguments:
      @ scheduler              <Scheduler> [Scheduler("sge", "qsub", "all.q")]
          : Job scheduler
      @ centromere_reads_fname <str>       ["centromere_reads.pkl"]
          : File of centromere reads
      @ out_fname              <str>       ["centromere_reads_unsync_overlaps.pkl"]
          : Output file name
      @ offset                 <int>       [1]
          : `offset` units around both boundaries of a read are not used as k-units.
      @ k_for_unit             <int>       [2]
          : `k_for_unit`-units will be used for initial mapping.
      @ k_for_spectrum         <int>       [13]
          : `k_for_spectrum`-mer spectrums of a `k_for_unit`-unit and a whole read sequence are
            initially compared.
      @ min_kmer_ovlp          <float>     [0.4]
          : Minimum ratio required for overlap of k-mer spectrums between a `k`-unit and a read.
      @ max_init_diff          <float>     [0.02]
          : Threshold for initial `k`-unit mapping.
      @ max_diff               <float>     [0.02]
          : Threshold for sequence dissimilarity of final overlaps.
    """
    n_distribute           : int
    n_core                 : int
    scheduler              : Scheduler = Scheduler("sge", "qsub", "all.q")
    centromere_reads_fname : str       = "centromere_reads.pkl"
    out_fname              : str       = "centromere_reads_unsync_overlaps.pkl"
    offset                 : int       = 1
    k_for_unit             : int       = 2
    k_for_spectrum         : int       = 13
    min_kmer_ovlp          : float     = 0.4
    max_init_diff          : float     = 0.02
    max_diff               : float     = 0.02

    def __post_init__(self):
        run_command(f"mkdir -p {out_dir}; rm -f {out_dir}/*")

    def run(self):
        jids = []
        for i in range(self.n_distribute):
            index = str(i + 1).zfill(int(np.log10(self.n_distribute) + 1))
            out_fname = f"{out_dir}/{out_prefix}.{index}.pkl"
            script_fname = f"{out_dir}/{scatter_prefix}.{index}.sh"

            script = ' '.join(map(str, ["python -m vca.overlapper.ava_unsync_reads",
                                        self.centromere_reads_fname,
                                        out_fname,
                                        self.n_distribute,
                                        self.n_core,
                                        self.offset,
                                        self.k_for_unit,
                                        self.k_for_spectrum,
                                        self.min_kmer_ovlp,
                                        self.max_init_diff,
                                        self.max_diff,
                                        i]))

            jids.append(self.scheduler.submit(script,
                                              script_fname,
                                              job_name="ava_unsync",
                                              log_fname=f"{log_prefix}.{index}",
                                              n_core=self.n_core))

        self.scheduler.submit("sleep 1s",
                              gather_fname,
                              job_name="ava_unsync_gather",
                              log_fname=log_prefix,
                              depend=jids,
                              wait=True)

        merged = []
        fnames = run_command(f"find {out_dir} -name '{out_prefix}.*' | sort").strip().split('\n')
        for fname in fnames:
            merged += load_pickle(fname)
        save_pickle(sorted(merged), self.out_fname)


def svs_overlap_mult(read_id_pairs,
                     offset, k_for_unit, min_kmer_ovlp, max_init_diff, max_diff):
    return [svs_overlap(reads[a_read_id], reads[b_read_id],
                        rc_reads[a_read_id], rc_reads[b_read_id],
                        offset, k_for_unit, min_kmer_ovlp, max_init_diff, max_diff,
                        read_forward_specs, read_boundary_specs)
            for a_read_id, b_read_id in read_id_pairs]


if __name__ == "__main__":
    """Only for internal usage."""
    p = argparse.ArgumentParser()
    p.add_argument("centromere_reads_fname", type=str)
    p.add_argument("out_fname", type=str)
    p.add_argument("n_distribute", type=int)
    p.add_argument("n_core", type=int)
    p.add_argument("offset", type=int)
    p.add_argument("k_for_unit", type=int)
    p.add_argument("k_for_spectrum", type=int)
    p.add_argument("min_kmer_ovlp", type=float)
    p.add_argument("max_init_diff", type=float)
    p.add_argument("max_diff", type=float)
    p.add_argument("index", type=int)
    args = p.parse_args()

    # Load all reads
    centromere_reads = load_pickle(args.centromere_reads_fname)
    centromere_reads_by_id = {read.id: read for read in centromere_reads}

    # List up read ID pairs assigned to this job
    read_id_pairs = [(a_read.id, b_read.id)
                     for a_read in centromere_reads
                     for b_read in centromere_reads
                     if a_read.id < b_read.id]
    unit_n = -(-len(read_id_pairs) // args.n_distribute)
    read_id_pairs = read_id_pairs[args.index * unit_n:(args.index + 1) * unit_n]

    # Precompute reads involved in this job and k-mer spectrums of the reads
    read_ids = set([read_id for read_id_pair in read_id_pairs for read_id in read_id_pair])
    global reads
    global rc_reads
    global read_forward_specs
    global read_boundary_specs
    reads = {read_id: centromere_reads_by_id[read_id] for read_id in read_ids}
    rc_reads = {read_id: revcomp_read(centromere_reads_by_id[read_id]) for read_id in read_ids}
    read_forward_specs = {read.id: seq_to_forward_kmer_spectrum(read.seq, k=args.k_for_spectrum)
                          for read in reads.values()}
    read_boundary_specs = {}
    for read in reads.values():
        for strand in (0, 1):
            if strand == 1:
                read = revcomp_read(read)
            # prefix boundary
            start, end = read.units[args.offset].start, read.units[args.offset + args.k_for_unit - 1].end
            read_boundary_specs[(read.id, strand, start, end)] = \
                seq_to_forward_kmer_spectrum(read.seq[start:end], k=args.k_for_spectrum)
            # suffix boundary
            start, end = read.units[-args.offset - args.k_for_unit].start, read.units[-args.offset - 1].end
            read_boundary_specs[(read.id, strand, start, end)] = \
                seq_to_forward_kmer_spectrum(read.seq[start:end], k=args.k_for_spectrum)

    # Divide into read_pairs for each core
    unit_n = -(-len(read_id_pairs) // args.n_core)
    read_id_pairs_list = [(read_id_pairs[i * unit_n:(i + 1) * unit_n],
                           args.offset, args.k_for_unit, args.min_kmer_ovlp,
                           args.max_init_diff, args.max_diff)
                          for i in range(args.n_core)]

    overlaps = set()
    with Pool(args.n_core) as pool:
        for ret_list in pool.starmap(svs_overlap_mult, read_id_pairs_list):
            for ret in ret_list:
                overlaps.update(ret)

    save_pickle(sorted(overlaps), args.out_fname)
