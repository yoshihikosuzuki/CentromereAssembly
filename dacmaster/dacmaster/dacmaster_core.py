import os
import random
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from interval import interval
from logzero import logger

from .dacmaster_io import load_unit_fasta

from BITS.utils import run_command, revcomp
from BITS.core import run_edlib


def take_consensus_cyclic(seqs, out_dir):
    """
    Return consensus sequence of the given sequences using Consed.
    First one will be the seed, and cyclic alignment is used in the mapping of other ones to it.
    This requires <out_dir> for a temporary place of the Consed input file.
    """

    if not os.path.isdir(out_dir):
        run_command(f"mkdir {out_dir}")
    out_fname = os.path.join(out_dir, "input.seq")   # temporary file for Consed input
    seq_writer = open(out_fname, 'w')

    # Output the first sequence as seed for consensus
    seed_seq = seqs[0]
    seq_writer.write(f"{seed_seq}\n")

    for seq in seqs[1:]:
        seq = seq * 2   # duplicated target sequence for a proxy of cyclic alignment
        alignment = run_edlib(seed_seq, seq, mode="glocal")
        #print(alignment["start"], alignment["end"], len(seq)/2)
        seq = seq[alignment["start"]:alignment["end"]]   # mapped area
        seq_writer.write(f"{seq}\n")

    seq_writer.close()
    return run_command(f"consed {out_fname}").replace('\n', '')


class Peak:
    def __init__(self, root_dir, index, N, unit_len, density, start_len, end_len, units):
        self.root_dir = root_dir
        self.index = index   # "<root_dir>/<fname_prefix>.<index>.fasta" is the result file
        self.N = N   # num of merged peaks inside this peak (1 if an independent peak)
        self.unit_len = unit_len   # unit length for each merged peak   # NOTE: list
        self.density = density   # density for each merged peak   # NOTE: list
        self.start_len = start_len   # min unit length in this peak
        self.end_len = end_len   # max unit length in this peak
        self.units = units   # longer than <start_len> and shorter than <end_len>   # NOTE: pandas dataframe

    def take_intra_consensus(self, n_units_threshold=5):
        """
        For each tandem repeat that has at least <n_units_threshold> units,
        this calculates intra-TR consensus sequence of the units.
        """

        self.units_consensus = {}   # intra-TR consensus units
        index = 0
        for read_id, df_read in self.units.groupby("read_id"):   # for each read
            for path_id, df_path in df_read.groupby("path_id"):   # for each TR
                if len(df_path) < n_units_threshold:   # only small number of units inside the TR
                    continue
                consensus_seq = take_consensus_cyclic(list(df_path["sequence"]), "tmp")
                #print(consensus_seq)
                if len(consensus_seq) == 0 or consensus_seq[0] == 'W' or consensus_seq[0] == '*':
                    # XXX: TODO: this is just a workaround for the output of consed "Warning: tile overlaps did not align well" or "** EXISTING". Fix it.
                    continue
                self.units_consensus[index] = [read_id,
                                               path_id,
                                               df_path["start"].iloc[0],   # TODO: maybe no need
                                               df_path["end"].iloc[-1],
                                               len(consensus_seq),
                                               consensus_seq]
                index += 1

        self.units_consensus = pd.DataFrame.from_dict(self.units_consensus,
                                                      orient="index",
                                                      columns=("read_id",
                                                               "path_id",
                                                               "start",
                                                               "end",
                                                               "length",
                                                               "sequence"))

    def propose_cluster(self, remaining_read_id, cluster_id, out_dir):
        """
        From the remaining consensus units, first randomly choose a unit and
        then collect all units within diameter of 10% sequence difference from
        it.
        """

        cluster_units = np.zeros(self.units_consensus.shpae[0], dtype=int)   # <cluster_id> + 1 for cluster members, otherwise 0
        read_id = random.choice(remaining_read_id)
        query = self.units_consensus["sequence"].iloc[read_id]
        
        if not os.path.isdir(out_dir):
            run_command(f"mkdir {out_dir}")
        out_fname = os.path.join(out_dir, "input.seq")   # temporary file for Consed input
        seq_writer = open(out_fname, 'w')
    
        # Output the first sequence as seed for consensus
        seq_writer.write(f"{query}\n")
    
        for j in remaining_read_id:
            target = self.units_consensus["sequence"].iloc[j]
            
            # forward
            seq_f = target * 2
            alignment_f = run_edlib(query, seq_f, mode="glocal")
    
            # reverse complement
            seq_rc = revcomp(seq_f)
            alignment_rc = run_edlib(query, seq_rc, mode="glocal")
    
            if alignment_f["diff"] <= alignment_rc["diff"]:
                alignment = alignment_f
                seq = seq_f
            else:
                alignment = alignment_rc
                seq = seq_rc
    
            if alignment["diff"] < 0.1:   # 0.05とかにしてもいいかも
                cluster_units[j] = cluster_id + 1
                seq = seq[alignment["start"]:alignment["end"]]   # mapped area
                seq_writer.write(f"{seq}\n")
    
        seq_writer.close()
        return (cluster_units,
                run_command(f"consed {out_fname}").replace('\n', ''))

    def fix_cluster_assignment(self, remaining_read_id, cluster_id, consensus_seq):   # TODO: can I merge this and above one into a single method?
        """
        Given a set of consensus units proposed, return an actual cluster
        which has diameter of 10% sequence difference from consensus of the
        consensus units. The cluster size might become very smaller than the
        proposed one.
        """

        cluster_units = np.zeros(N, dtype=int)   # 1 for cluster members, otherwise 0
        read_id = random.choice(remaining_read_id)
        query = self.units_consensus["sequence"].iloc[read_id]
    
        for j in remaining_read_id:
            target = self.units_consensus["sequence"].iloc[j]
    
            # forward
            seq_f = target * 2
            alignment_f = run_edlib(query, seq_f, mode="glocal")
    
            if alignment_f["diff"] < 0.1:
                cluster_units[j] = cluster_id + 1   # add 1 to neutralize the default value -1 in the <assignment>
            else:
                # reverse complement
                seq_rc = revcomp(seq_f)
                alignment_rc = run_edlib(query, seq_rc, mode="glocal")
                if alignment_rc["diff"] < 0.1:
                    cluster_units[j] = cluster_id + 1
                    
        return cluster_units

    def find_representatives(self):
        """
        From all units in a peak, find representative monomers that do not align well to each other.
        This task includes start-position adjustment and phase-adjustment as well.

        In this step, determination of seed units is essential.
        """

        # Calculate consensus of the units for each TR
        self.take_intra_consensus()

        # Clustering of the consensus untis
        # Since what we are doing for now is just determine rough representative monomers,
        # this task is rather just to eliminate redundant consensus units.
        self.assignment = np.full(self.units_consensus.shape[0], dtype=int)   # cluster id assignment for each consensus unit
        self.clusters = {}   # mater unit sequences
        cluster_id = 0
        remaining_read_id = np.arange(N)   # array of indices of consensus units still remaining

        while True:
            cluster_units, consensus_seq = self.propose_cluster(remaining_read_id,
                                                                cluster_id,
                                                                f"{self.root_dir}tmp/")
            
            # XXX: consensus sequence may be "EXITING", "** WARNING", or Nothing!!!!
            # TODO: maybe I should first debug Consed
            
            cluster_size = cluster_units[cluster_units > 0].shape[0]
            print(cluster_id, cluster_size)
            
            #print(consensus_seq)
            if cluster_size >= self.units_consensus.shape[0] * 0.01:
                print("accepted")
                
                # Re-align the consensus sequence to the remaining consensus units
                cluster_units = self.fix_cluster_assignment(remaining_read_id,
                                                            cluster_id,
                                                            consensus_seq)   # TODO: more efficient way?
                
                cluster_size = cluster_units[cluster_units > 0].shape[0]
                print("actual cluster size =", cluster_size)
                
                if cluster_size >= self.units_consensus.shape[0] * 0.01:   # filter by cluster size again
                    print("accepted")
                    self.clusters[cluster_id] = consensus_seq
                    self.assignment += cluster_units
                    cluster_id += 1
                else:
                    print("rejected")
            else:
                print("rejected")
                self.assignment -= cluster_units
                
            remaining_read_id = np.where(self.assignment == -1)[0]
            print("# of remaining units =", remaining_read_id.shape[0])
            #print(clusters)
            if remaining_read_id.shape[0] < self.units_consensus.shape[0] * 0.1:
                break
        
        self.assignment[np.where(self.assignment < 0)] = -1

        logger.info(f"{len(self.clusters)} representative unit candidates:\n{self.clusters}")
        self.clusters = pd.DataFrame.from_dict(self.clusters, orient="index")
        # TODO: add cluster size column and sort by it

        # Remove redundant untis which are similar to another one
        # After that, the remaining centroids are representative units
        dm = np.zeros((self.clusters.shape[0], self.clusters.shape[0]), dtype=float)
        for i in range(self.clusters.shape[0] - 1):
            for j in range(1, self.clusters.shape[0]):
                #print(f"# {i} vs {j} (f)")
                alignment = run_edlib(self.clusters[0][i], self.clusters[0][j] * 2, mode='glocal')
                f = alignment["diff"]
                #print(f"# {i} vs {j} (rc)")
                alignment = run_edlib(self.clusters[0][i], revcomp(self.clusters[0][j] * 2), mode='glocal')
                rc = alignment["diff"]
                
                dm[i, j] = dm[j, i] = min([f, rc])

        self.representative_units = {self.clusters[0][0]}
        for i in range(1, self.clusters.shape[0]):
            flag_add = 1
            for j in range(i):
                if dm[i, j] < 0.1:
                    flag_add = 0
                    break
            if flag_add == 1:
                self.representative_units.add(clusters[0][i])

        logger.info(f"{len(self.representative_units)} after redundancy removal:\n{self.representative_units}")

        for i, r_unit in enumerate(self.representative_units):
            print(f">representative/{i}/0_{len(r_unit)}\n{r_unit}\n")   # add cluster_size information
        

class Runner:
    def __init__(self,
                 unit_fasta,
                 peaks_dir,
                 min_len=50,   # peak length must be longer than this
                 max_len=1000,   # must be shorter than this
                 band_width=5,   # parameter for KDE
                 min_density=0.001,   # threshold of peak hight
                 deviation=0.1):   # units inside of "peak_len +- deviation %" are collected as Peak.units

        self.units = load_unit_fasta(unit_fasta)
        self.peaks_dir = peaks_dir
        self.min_len = min_len
        self.max_len = max_len
        self.band_width = band_width
        self.min_density = min_density
        self.deviation = deviation

    def run(self):
        self.smooth_unit_len_dist()
        self.detect_peaks()

    def smooth_unit_len_dist(self):
        """
        Calculate unit length distribution smoothed by kernel density estimation.
        """
        
        self.unit_lens = np.array(self.units["length"])

        self.ul = self.unit_lens[(self.min_len < self.unit_lens) & (self.unit_lens < self.max_len)]
        self.ls = np.linspace(self.min_len, self.max_len, self.max_len - self.min_len + 1, dtype=int)

        self.kde = KernelDensity(kernel='gaussian', bandwidth=self.band_width).fit(self.ul.reshape(-1, 1))
        self.dens = np.exp(self.kde.score_samples(self.ls.reshape(-1, 1)))

    def detect_peaks(self):
        """
        Detect peaks in the unit length distribution by simple sweep line.
        Adjascent peaks close to each other will be automatically merged.
        """

        self.peak_intervals = interval()
        self.peak_info = []
        prev_n_peaks = 0
        for i in range(1, len(self.dens) - 1):
            if ((self.dens[i] >= self.min_density) and
                (self.dens[i - 1] < self.dens[i]) and
                (self.dens[i] > self.dens[i + 1])):

                self.peak_intervals |= interval[-(- self.ls[i] * (1. - self.deviation) // 1),
                                           int(self.ls[i] * (1. + self.deviation))]
                peak_info = (self.ls[i], self.dens[i])

                logger.info(f"Peak detected: length = {self.ls[i]} bp, density = {self.dens[i]}")
                if prev_n_peaks == len(self.peak_intervals):
                    self.peak_info[-1].append(peak_info)
                    logger.info("Merged to the previous peak.")
                else:   # new peak interval is independent
                    self.peak_info.append([peak_info])
                    prev_n_peaks += 1
    
        # For each peak interval, create Peak class instance
        self.peaks = []
        for i, intvl in enumerate(self.peak_intervals.components):
            peak_info = self.peak_info[i]
            peak_n = len(peak_info)   # num of peaks merged
            peak_len, peak_dens = list(zip(*peak_info))   # list for each merged peak
            peak_start, peak_end = intvl[0]   # min peak_len - deviation, max peak_len + deviation
            peak_units = (self.units[self.units["length"] >= peak_start]
                          .pipe(lambda df: df[df["length"] <= peak_end]))
            self.peaks.append(Peak(self.peaks_dir,
                                   i,
                                   peak_n,
                                   peak_len,
                                   peak_dens,
                                   peak_start,
                                   peak_end,
                                   peak_units))
