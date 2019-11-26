from .core import gen_unique_seq, mutate_seq


def is_tandem(seq):
    """Check if `seq` is tandemly repetitive or not, assuming it is error-free."""
    L = len(seq)
    for i in range(1, int(L / 2) + 1):
        if L % i != 0:
            continue
        if seq == seq[:i] * int(L / i):
            return True
    return False


def gen_random_non_tandem_repeat_seq(length):
    """Generate a NON tandem repeat sequence of `length` bp, which can be a unit sequence of
    a tandem repeat."""
    while True:
        seq = gen_unique_seq(length)
        if not is_tandem(seq):
            return seq


def gen_random_tandem_repeats(n_seq, unit_len, n_copy, edit_weights, prefix_len, suffix_len):
    """Generate `n_seq` sequences where each consists of:

      [unique sequence of `prefix_len` bp]
      + [(unique sequence of `unit_len` bp) * `n_copy`]
      + [unique sequence of `suffix_len` bp]

    `edit_weights` is a list of weights for each edit operation, {=, X, I, D} in this order.
    Ex.) (85, 1, 11, 3) for PacBio CLR
    """
    def gen_tandem_repeat():
        pre_seq = gen_random_non_tandem_repeat_seq(prefix_len)
        suf_seq = gen_random_non_tandem_repeat_seq(suffix_len)
        unit_seq = gen_random_non_tandem_repeat_seq(unit_len)
        tr_seq = unit_seq * n_copy
        entire_seq = pre_seq + tr_seq + suf_seq
        return (mutate_seq(entire_seq, edit_weights), unit_seq)

    return [gen_tandem_repeat() for i in range(n_seq)]
