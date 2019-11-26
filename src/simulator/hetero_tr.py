from .core import gen_unique_seq, mutate_seq, sequence_seq


def gen_hetero_tandem_repeat(unit_length, copy_num, mutate_rate=1, flanking_length=10000):
    unit = gen_unique_seq(unit_length)
    seq = unit
    for i in range(copy_num - 1):
        unit = sequence_seq(unit, (100 - mutate_rate, mutate_rate / 3, mutate_rate / 3, mutate_rate / 3))
        seq += unit
    return gen_unique_seq(flanking_length) + seq + gen_unique_seq(flanking_length)
