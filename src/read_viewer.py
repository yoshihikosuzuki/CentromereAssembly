from copy import deepcopy
from dataclasses import dataclass
from BITS.seq.plot import DotPlot
from BITS.plot.plotly import make_line, make_scatter, make_layout, show_plot
from BITS.util.proc import run_command
from .datruf.io import load_dumps, find_inner_alignments
from .types import TRRead


@dataclass(eq=False)
class ReadViewer:
    db_fname  : str
    las_fname : str
    gepard    : str
    out_dir   : str = "tmp"
    plot_size : int = 1000

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")

    def show(self, in_read, dot_plot=False, alignment_plot=True, out_html=None):
        """<read> must be ID of a read in <db_fname> or TRRead object instance."""
        assert isinstance(in_read, (int, TRRead)), "Unsupported read type"
        read = deepcopy(in_read) if isinstance(TRRead) else TRRead(id=in_read)
        assert read.id is not None, "Read ID must be set"

        if dot_plot:
            self._dot_plot(read)
        if alignment_plot:
            self._alignment_plot(read, out_html)

    def _dot_plot(self, read):
        out_fasta = f"{self.out_dir}/{read.id}.fasta"
        run_command(f"DBshow {self.db_file} {read.id} > {out_fasta}")
        DotPlot(self.gepard, self.out_dir).plot_fasta(out_fasta, out_fasta)

    def _alignment_plot(self, read, out_html):
        read_len = calc_read_len(read.id, self.db_fname)   # read length

        # Load information of TR intervals and self alignments
        # and create a new TRRead object having enough data to make a plot
        dump = load_dumps(read.id, read.id, self.db_fname, self.las_fname)[0]
        inner_alignments = find_inner_alignments(dump)

        # Create line objects which represent the read itself, TR intervals, and self alignments
        shapes = [make_line(0, 0, read_len, read_len, "grey", 3)]   # read
        shapes += [make_line(tr.start, tr.start, tr.end, tr.end, "black", 3) for tr in dump.trs]   # TRs
        shapes += [make_line(aln.ab, aln.bb, aln.ae, aln.be,   # self alignments
                             col=("purple" if aln in inner_alignments
                                  else "black" if 0.95 <= aln.slope <= 1.05
                                  else "yellow"),   # abnormal slope (= noisy)
                             width=3 if aln in inner_alignments else 1)
                   for aln in dump.alignments]
            
        # Create traces of start/end positions of the self alignments and TR intervals
        trace_start = make_scatter([aln.ab for aln in dump.alignments],
                                   [aln.bb for aln in dump.alignments],
                                   name="start")
        trace_end   = make_scatter([aln.ae for aln in dump.alignments],
                                   [aln.be for aln in dump.alignments],
                                   name="end")
        trace_tr    = make_scatter([tr.start for tr in dump.trs],
                                   [tr.start for tr in dump.trs],
                                   text=list(range(1, len(dump.trs) + 1)),
                                   text_pos="top right", text_size=10, text_col="grey", mode="text",
                                   name="TR interval")
        traces = [trace_start, trace_end, trace_tr]

        # Add shapes and traces if <read> has data on units
        if read.units is not None:
            shapes += [make_line(unit.start, unit.start, unit.end, unit.end, "black", 5)
                       for unit in read.units]   # on diagonal   # TODO: change color based on the attributes
            trace_unit = make_scatter([unit.start for unit in read.units],
                                      [unit.start for unit in read.units],
                                      text=[f"{i} " for i in range(len(read.units))],
                                      text_pos="bottom left", text_size=10, text_col="black", mode="text",
                                      name="TR unit")
            trace_info = make_scatter([unit.start for unit in read.units],
                                      [unit.start for unit in read.units],
                                      text=[f"unit {i}<br>{unit.start}:{unit.end}"
                                            for i, unit in enumerate(read.units)],
                                      col="black",
                                      show_legend=False)
            traces += [trace_unit, trace_info]

        # Show plot
        layout = make_layout(self.plot_size, self.plot_size, title=f"{read.id}",
                             x_range=[-read_len * 0.05, read_len + 100],
                             y_range=[0, read_len],
                             x_grid=False, y_grid=False, y_reversed=True, shapes=shapes)
        show_plot(traces, layout, out_html)


def calc_read_len(read_id, db_fname):
    """Call DBdump with <read_id> and extract read length information from the output."""
    command = f"DBdump -s {db_fname} {read_id}"
    return int(run_command(command).strip().split('\n')[-1].split(' ')[1])
