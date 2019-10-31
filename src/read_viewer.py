from dataclasses import dataclass
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex, XKCD_COLORS
from logzero import logger
from BITS.clustering.seq import ClusteringSeq
from BITS.seq.io import save_fasta
from BITS.seq.align import EdlibRunner
from BITS.seq.dot_plot import DotPlot
from BITS.plot.plotly import make_line, make_rect, make_scatter, make_layout, show_plot
from BITS.util.proc import run_command
from .datruf.io import load_tr_reads
from .datruf.find_units import find_inner_alignments
from .types import TRRead


@dataclass(eq=False)
class TRReadViewer:
    """Class for plotting reads with tandem repeats in Jupyter Notebook.

    Basic usage (in Jupyter Notebook):
      > from vca import TRReadViewer
      > v = ReadViewer(db_fname, las_fname)
      > v.show(read)


    positional arguments of `TRReadViewer`:
      @ db_fname  <str> : DAZZ_DB `.db` file of the reads.
      @ las_fname <str> : `TAN.*.las` 

    optional arguments of `TRReadViewer`:
      @ gepard  <str> [None]  : Command to run Gepard. It would be like `f"java -cp {gepard_jar}
                                org.gepard.client.cmdline.CommandLine -matrix {gepard_mat}"`.
                                This is required if `dot_plot=True` in the `show()` method.
      @ out_dir <str> ["tmp"] : Directory for outputting intermediate files.


    positional arguments of the `show()` method:
      @ a_read <int|TRRead> : Read ID used in `db_fname` (int) or a `vca.types.TRRead` object.

    optional arguments of the `show()` method:
      @ b_read         <int|TRRead> [None]
          : If specified, show not self-vs-self plots but `a_read`-vs-`b_read` plots.
      @ dot_plot       <bool>       [False]
          : If `True`, show a dot plot with Gepard. `gepard` must not be `None`.
      @ alignment_plot <bool>       [True]
          : If `True`, show an alignment plot.
      @ plot_size      <int>        [700]
          : Size of the alignment plot.
      @ out_fname      <str>        [None]
          : If specified, alignment plot is output to this file.
    """
    db_fname  : str
    las_fname : str
    gepard    : str = None
    out_dir   : str = "tmp"

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")

    def show(self, a_read, b_read=None, dot_plot=False, alignment_plot=True, plot_size=700, out_fname=None):
        a_read = self.load_read(a_read)
        if b_read is not None:
            b_read = self.load_read(b_read)

        if dot_plot:
            self._dot_plot(a_read, b_read)
        if alignment_plot:
            if b_read is None:
                self._alignment_plot_self(a_read, plot_size, out_fname)
            else:
                self._alignment_plot_other(a_read, b_read, plot_size, out_fname)

    def load_read(self, read):
        assert isinstance(read, (int, TRRead)), "Type of `read` must be int or TRRead"

        if isinstance(read, int):   # read ID is specified
            read = load_tr_reads(read, read, self.db_fname, self.las_fname, return_all=True)[0]
        assert read.id is not None, "`read.id` must be specified"
        return read

    def read_to_fasta(self, read):
        out_fasta = f"{self.out_dir}/{read.id}.fasta"
        save_fasta({f"{read.id}{'' if read.strand == 0 else '_rc'}": read.seq},
                   out_fasta, sort=False, width=100)
        return out_fasta

    def _dot_plot(self, a_read, b_read):
        assert self.gepard is not None, "`gepard` must be specified"

        a_out_fasta = self.read_to_fasta(a_read)
        if b_read is None:   # self-vs-self of `a_read`
            DotPlot(self.gepard, self.out_dir).plot_fasta(a_out_fasta, a_out_fasta)
        else:   # `a_read` vs `b_read`
            b_out_fasta = self.read_to_fasta(b_read)
            DotPlot(self.gepard, self.out_dir).plot_fasta(a_out_fasta, b_out_fasta)

    def _alignment_plot_self(self, read, plot_size, out_fname):
        traces, shapes = [], []

        # Read as a shape on the diagonal
        shapes += [make_line(0, 0, read.length, read.length, "grey", 3)]

        # Add shapes and traces of TR intervals and self alignments if exist
        if read.trs is not None and read.alignments is not None:
            inner_alignments = find_inner_alignments(read)

            shapes += [make_line(tr.start, tr.start, tr.end, tr.end, "black", 3)
                       for tr in read.trs]   # TRs
            shapes += [make_line(aln.ab, aln.bb, aln.ae, aln.be,   # self alignments
                                 col=("purple" if aln in inner_alignments
                                      else "black" if 0.95 <= aln.slope <= 1.05
                                      else "yellow"),   # abnormal slope (= noisy)
                                 width=3 if aln in inner_alignments else 1)
                       for aln in read.alignments]

            # Create traces of start/end positions of the self alignments and TR intervals
            trace_start = make_scatter([aln.ab for aln in read.alignments],
                                       [aln.bb for aln in read.alignments],
                                       name="start")
            trace_end   = make_scatter([aln.ae for aln in read.alignments],
                                       [aln.be for aln in read.alignments],
                                       name="end")
            trace_tr    = make_scatter([tr.start for tr in read.trs],
                                       [tr.start for tr in read.trs],
                                       text=list(range(1, len(read.trs) + 1)),
                                       text_pos="top right", text_size=10, text_col="grey", mode="text",
                                       name="TR interval")
            traces += [trace_start, trace_end, trace_tr]

        # Add shapes and traces of TR units if exist
        if read.units is not None:
            logger.warn(f"Distances are inaccurate for unsynchronized units.")
            d_to_c = np.vectorize(lambda x: rgb2hex(cm.Blues_r(x)))

            # Sequence dissimilarity between the raw units
            c_raw = ClusteringSeq(read.unit_seqs, revcomp=False, cyclic=False if read.synchronized else True)
            c_raw.calc_dist_mat()
            raw_dist = [[round(c_raw.s_dist_mat[i][j] * 100, 2) for j in range(len(read.units))]
                        for i in range(len(read.units))]

            # Sequence dissimilarity between repr units assigned to the units
            if read.synchronized:   # units must be encoded with repr units if synchronized
                c_repr = ClusteringSeq([read.repr_units[unit.repr_id] for unit in read.units],
                                       revcomp=False, cyclic=False)
                c_repr.calc_dist_mat()
                repr_dist = [[round(c_repr.s_dist_mat[i][j] * 100, 2) for j in range(len(read.units))]
                             for i in range(len(read.units))]
                cols = d_to_c(c_repr.s_dist_mat / np.max(c_repr.s_dist_mat))
            else:
                repr_dist = [['-' for j in range(len(read.units))] for i in range(len(read.units))]
                cols = d_to_c(c_raw.s_dist_mat / np.max(c_raw.s_dist_mat))

            repr_ids = [unit.repr_id if unit.repr_id is not None else '-' for unit in read.units]
            
            unit_texts = [f"Unit {i} (repr={repr_ids[i]}) vs Unit {j} (repr={repr_ids[j]})<br>"
                          f"{raw_dist[i][j]}% diff (raw), {repr_dist[i][j]}% diff (repr)"
                          for i in range(len(read.units))
                          for j in range(i + 1, len(read.units))]
            unit_cols = [cols[i][j]
                         for i in range(len(read.units))
                         for j in range(i + 1, len(read.units))]

            # Distance matrix
            shapes += [make_rect(read.units[i].start, read.units[j].start,
                                 read.units[i].end, read.units[j].end,
                                 fill_col=cols[i][j])
                       for i in range(len(read.units)) for j in range(i + 1, len(read.units))]

            traces += [make_scatter([((read.units[i].start + read.units[i].end) / 2)
                                     for i in range(len(read.units))
                                     for j in range(i + 1, len(read.units))],
                                    [((read.units[j].start + read.units[j].end) / 2)
                                     for i in range(len(read.units))
                                     for j in range(i + 1, len(read.units))],
                                    text=unit_texts, col=unit_cols,
                                    marker_size=3, show_scale=True, show_legend=False)]
            
            # Units on diagonal
            id_to_col = {}
            for i, col in enumerate(XKCD_COLORS.values()):   # NOTE: up to ~1000 cols
                id_to_col[i] = col
            shapes += [make_line(unit.start, unit.start, unit.end, unit.end,
                                 id_to_col[unit.repr_id] if read.synchronized else "black",
                                 5)
                       for unit in read.units]
            trace_unit = make_scatter([unit.start for unit in read.units],
                                      [unit.start for unit in read.units],
                                      text=[f"{i} " for i in range(len(read.units))],
                                      text_pos="bottom left", text_size=10, text_col="black",
                                      mode="text", name="TR units")

            er_global = EdlibRunner("global", revcomp=False, cyclic=False)
            diff_from_repr = [round(100 * er_global.align(read.repr_units[unit.repr_id],
                                                          read.seq[unit.start:unit.end]).diff, 2)
                              if read.synchronized else '-'
                              for i, unit in enumerate(read.units)]

            texts = [f"Unit {i} (repr={repr_ids[i]}; strand={unit.strand if read.synchronized else '-'})<br>"
                     f"[{unit.start}:{unit.end}] ({unit.length} bp)<br>"
                     f"{diff_from_repr[i]}% diff from repr unit"
                     for i, unit in enumerate(read.units)]
            
            trace_info = make_scatter([unit.start for unit in read.units],
                                      [unit.start for unit in read.units],
                                      text=texts,
                                      col=([id_to_col[unit.repr_id] for unit in read.units]
                                           if read.synchronized else "black"),
                                      show_legend=False)
            traces += [trace_unit, trace_info]

        # Show plot
        layout = make_layout(plot_size * 1.05, plot_size, title=f"Read {read.id} (strand={read.strand})",
                             x_range=[-read.length * 0.05, read.length + 100],
                             y_range=[0, read.length],
                             x_grid=False, y_grid=False, y_reversed=True, shapes=shapes)
        layout.update(dict(margin=dict(l=0, t=30, b=0)))
        show_plot(traces, layout, out_fname)

    def _alignment_plot_other(self, a_read, b_read, plot_size, out_fname):
        assert a_read.units is not None and b_read.units is not None, "`[a|b]_read.units` must not be None"
        assert a_read.synchronized and b_read.synchronized, "Both reads must be synchronized"

        id_to_col = {}   # {repr_id: color}
        for i, col in enumerate(XKCD_COLORS.values()):   # NOTE: up to ~1000 cols
            id_to_col[i] = col

        traces, shapes = [], []

        # Unit encodings
        er_global = EdlibRunner("global", revcomp=False, cyclic=False)
        a_diff_from_repr = [round(100 * er_global.align(a_read.repr_units[unit.repr_id],
                                                        a_read.seq[unit.start:unit.end]).diff, 2)
                            for unit in a_read.units]
        b_diff_from_repr = [round(100 * er_global.align(b_read.repr_units[unit.repr_id],
                                                        b_read.seq[unit.start:unit.end]).diff, 2)
                            for unit in b_read.units]

        shapes += [make_line(0, -b_read.length * 0.01, a_read.length, -b_read.length * 0.01, "grey", 3),
                   make_line(-a_read.length * 0.01, 0, -a_read.length * 0.01, b_read.length, "grey", 3)]
        shapes += [make_line(unit.start, -b_read.length * 0.01, unit.end, -b_read.length * 0.01,
                             id_to_col[unit.repr_id], 5)
                   for unit in a_read.units]
        shapes += [make_line(-a_read.length * 0.01, unit.start, -a_read.length * 0.01, unit.end,
                             id_to_col[unit.repr_id], 5)
                   for unit in b_read.units]
        traces += [make_scatter([(unit.start + unit.end) / 2 for unit in a_read.units],
                                [-b_read.length * 0.01 for unit in a_read.units],
                                text=[f"Unit {i} (repr={unit.repr_id})"
                                      f"[{unit.start}:{unit.end}] ({unit.length} bp)<br>"
                                      f"{a_diff_from_repr[i]}% diff from repr unit"
                                      for i, unit in enumerate(a_read.units)],
                                col=[id_to_col[unit.repr_id] for unit in a_read.units],
                                show_legend=False),
                   make_scatter([-a_read.length * 0.01 for unit in b_read.units],
                                [(unit.start + unit.end) / 2 for unit in b_read.units],
                                text=[f"Unit {i} (repr={unit.repr_id})"
                                      f"[{unit.start}:{unit.end}] ({unit.length} bp)<br>"
                                      f"{b_diff_from_repr[i]}% diff from repr unit"
                                      for i, unit in enumerate(b_read.units)],
                                col=[id_to_col[unit.repr_id] for unit in b_read.units],
                                show_legend=False)]

        # Distance matrix
        er_global = EdlibRunner("global", revcomp=False, cyclic=False)
        raw_dist = np.array([[er_global.align(a_unit_seq, b_unit_seq).diff
                              for b_unit_seq in b_read.unit_seqs]
                             for a_unit_seq in a_read.unit_seqs],
                            dtype=np.float32)
        repr_dist = np.array([[er_global.align(a_read.repr_units[a_unit.repr_id],
                                               b_read.repr_units[b_unit.repr_id]).diff
                               for b_unit in b_read.units]
                              for a_unit in a_read.units],
                             dtype=np.float32)
        d_to_c = np.vectorize(lambda x: rgb2hex(cm.Blues_r(x)))
        cols = d_to_c(repr_dist / np.max(repr_dist))

        # Heatmap of the distance matrix
        text = [f"Read {a_read.id} unit {i} (repr={a_read.units[i].repr_id}) vs "
                f"Read {b_read.id} unit {j} (repr={b_read.units[j].repr_id})<br>"
                f"{round(100 * raw_dist[i][j], 2)}% diff (raw), "
                f"{round(100 * repr_dist[i][j], 2)}% diff (repr)"
                for i in range(len(a_read.units)) for j in range(len(b_read.units))]
        col = [cols[i][j] for i in range(len(a_read.units)) for j in range(len(b_read.units))]

        shapes += [make_rect(a_read.units[i].start, b_read.units[j].start,
                             a_read.units[i].end, b_read.units[j].end,
                             fill_col=cols[i][j])
                   for i in range(len(a_read.units)) for j in range(len(b_read.units))]

        traces += [make_scatter([(a_read.units[i].start + a_read.units[i].end) / 2
                                 for i in range(len(a_read.units))
                                 for j in range(len(b_read.units))],
                                 [(b_read.units[j].start + b_read.units[j].end) / 2
                                  for i in range(len(a_read.units))
                                  for j in range(len(b_read.units))],
                                  text=text, col=col, marker_size=3,
                                show_scale=True, show_legend=False)]

        layout = make_layout(plot_size * 1.05, plot_size,
                             x_title=f"Read {a_read.id} (strand={a_read.strand})",
                             y_title=f"Read {b_read.id} (strand={b_read.strand})",
                             x_range=(-a_read.length * 0.05, a_read.length),
                             y_range=(0, b_read.length),
                             x_grid=False, y_grid=False, x_zeroline=False, y_zeroline=False,
                             y_reversed=True, shapes=shapes)
        layout["yaxis"]["scaleanchor"] = "x"
        layout["margin"] = dict(l=50, t=0, b=50)

        show_plot(traces, layout, out_fname)
