from dataclasses import dataclass
from BITS.seq.plot import DotPlot
from BITS.plot.plotly import make_line, make_scatter, make_layout, show_plot
from BITS.util.proc import run_command
from .datruf import load_dumps, filter_alignments


@dataclass(repr=False, eq=False)
class ReadViewer:
    db_file: str
    las_file: str
    out_dir: str
    gepard: str

    def __post_init__(self):
        run_command(f'mkdir -p {self.out_dir}')

    def show(self, read, dot_plot=False, alignment_plot=True, size=1000, out_html=None):
        if dot_plot:
            self._dot_plot(read)
        if alignment_plot:
            self._alignment_plot(read, size, out_html)

    def _alignment_plot(self, read, size=1000, out_html=None):
        # Load TR intervals and self alignments
        tr_intervals, alignments = load_dumps(read.id, read.id, self.db_file, self.las_file)
        tr_intervals, alignments = tr_intervals[read.id], alignments[read.id]
        tr_alignments = filter_alignments(tr_intervals, alignments)

        read_len = len(read.seq)
        shapes = [make_line(0, 0, read_len, read_len, 'grey', 3)]   # diagonal (= read)
        for start, end in tr_intervals:   # on diagonal
            shapes.append(make_line(start, start, end, end, 'black', 3))

        for ab, ae, bb, be, distance, slope in alignments:
            if (ab, ae, bb, be) in tr_alignments:
                col, width = 'purple', 3
            else:
                col = 'black' if 0.95 <= slope <= 1.05 else 'yellow'
                width = 1
            shapes.append(make_line(ab, bb, ae, be, col, width))   # self alignment

        traces = [make_scatter([x[0] for x in alignments],   # ab
                               [x[2] for x in alignments],   # bb
                               name='start'),
                  make_scatter([x[1] for x in alignments],   # ae
                               [x[3] for x in alignments],   # be
                               name='end'),
                  make_scatter([x[0] for x in tr_intervals],
                               [x[0] for x in tr_intervals],
                               text=list(range(1, len(tr_intervals) + 1)),
                               text_pos='top right', text_size=10, text_col='grey', mode='text',
                               name='TR interval')]

        if hasattr(read, 'units'):
            shapes += [make_line(unit.start, unit.start, unit.end, unit.end, 'black', 5)
                       for unit in read.units]
            traces += [make_scatter([unit.start for unit in read.units],
                                    [unit.start for unit in read.units],
                                    text=[f'{i} ' for i in range(len(list(read.units)))],
                                    text_pos='bottom left', text_size=10, text_col='black', mode='text',
                                    name='TR unit'),
                       make_scatter([unit.start for unit in read.units],
                                    [unit.start for unit in read.units],
                                    text=[f'unit {i}<br>{unit.start}:{unit.end}'
                                          for i, unit in enumerate(read.units)],
                                    col='black',
                                    show_legend=False)]

        layout = make_layout(size, size, title=f'{read.id}',
                             x_range=[-read_len * 0.05, read_len + 100],
                             y_range=[0, read_len],
                             x_grid=False, y_grid=False, y_reversed=True, shapes=shapes)
        
        show_plot(traces, layout, out_html)

    def _dot_plot(self, read):
        out_fasta = f'{self.out_dir}/{read.id}.fasta'
        run_command(f'DBshow {self.db_file} {read.id} > {out_fasta}')
        DotPlot(self.gepard, self.out_dir).plot_fasta(out_fasta, out_fasta)
