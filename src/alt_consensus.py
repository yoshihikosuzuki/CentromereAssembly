from collections import Counter, defaultdict
from BITS.seq.align import EdlibRunner


class PairwiseAlignment:
    def __init__(self, a_seq, b_seq):
        er = EdlibRunner("global", revcomp=False, cyclic=False)
        self.fcigar = er.align(b_seq.lower(), a_seq.lower()).cigar.flatten().string   # NOTE: b vs a; be careful!
        self.source, self.target = '', ''
        s_pos, t_pos = 0, 0
        for c in self.fcigar:
            if c == '=' or c == 'X':
                self.source += a_seq[s_pos]
                self.target += b_seq[t_pos]
                s_pos += 1
                t_pos += 1
            elif c == 'I':
                self.source += '-'
                self.target += b_seq[t_pos]
                t_pos += 1
            else:
                self.source += a_seq[s_pos]
                self.target += '-'
                s_pos += 1
        
    def show(self, by_cigar=False):
        if by_cigar:   # standard alignment like BLAST
            print(self.source)
            print(self.fcigar)
            print(self.target)
        else:
            print(''.join([' ' if c == '=' else self.source[i] for i, c in enumerate(self.fcigar)]))
            print(''.join([self.source[i] if c == '=' else ' ' for i, c in enumerate(self.fcigar)]))
            print(''.join([' ' if c == '=' else self.target[i] for i, c in enumerate(self.fcigar)]))


def count_variants(cluster_cons_unit, cluster_units):
    """Given a set of unit sequences <units> in a cluster, calculate the composition of
    nucleotides including '-' (= distribution of each )
    for each position on <cluster_cons_unit> as a seed.
    from which <units> are generated, compute the variations (= nucleotides inconsistent between
    <units> and <cluster_cons_unit> and their relative frequency).
    Since a cluster should be homogeneous (i.e., mono-source), the relative frequencies are
    expected to be not much larger than sequencing error.
    """
    assert cluster_cons_unit != "", "Empty strings are not allowed"
    # TODO: how to decide "same variant?" especially for multiple variations on same position (but slightly different among units)?
    variants = Counter()
    for unit in cluster_units:
        assert unit != "", "Empty strings are not allowed"
        alignment = PairwiseAlignment(cluster_cons_unit, unit)   # alignment.fcigar(cluster_cons_unit) = unit
        tpos = 0
        var_index = 0   # positive values for continuous insertions
        for i, c in enumerate(alignment.fcigar):
            if c == '=':
                var_index = 0
            elif c == 'I':
                var_index += 1
            if c != '=':
                variants[(tpos, var_index, c, alignment.target[i])] += 1   # TODO: multiple D on the same pos are aggregated
            if c != 'I':
                tpos += 1
        assert tpos == len(cluster_cons_unit)
    return variants


def consensus_alt(in_seqs, seed_choice="original"):
    """Compute a consensus sequence among `seqs: List[str]` by a simple majority vote for each position
    of the alignment pileup that is made by globally aligning a seed sequence and each of the other sequences.
    """
    
    # Choose the seed and move it to the first element of `in_seqs`
    assert seed_choice in ("original", "median", "longest"), "Invalid `seed_choise`"
    if seed_choice != "original":
        in_seqs = sorted(in_seqs, key=lambda x: len(x), reverse=True)   # "longest"
        if seed_choice == "median":
            index = len(in_seqs) // 2
            in_seqs = [in_seqs[index]] + in_seqs[:index] + in_seqs[index + 1:]
            
    var_counts = count_variants(in_seqs[0], in_seqs[1:])

    freqs = defaultdict(list)
    for (pos, subpos, _type, base), count in var_counts.items():
        freqs[(pos, subpos)].append((base, count))
        
    cons = ""

    for pos in range(len(in_seqs[0]) + 1):
        subpos = 1
        # insertions
        while (pos, subpos) in freqs:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in freqs[(pos, subpos)]:
                base_counts[base] = count
            base_counts['-'] = len(in_seqs) - sum(base_counts.values())
            cons += base_counts.most_common()[0][0]
            subpos += 1
        if pos == len(in_seqs[0]):
            break
        # others
        if (pos, 0) not in freqs:
            cons += in_seqs[0][pos]
        else:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in freqs[(pos, 0)]:
                base_counts[base] = count
            base_counts[in_seqs[0][pos]] = len(in_seqs) - sum(base_counts.values())
            cons += base_counts.most_common()[0][0]

    return cons.replace('-', '')
