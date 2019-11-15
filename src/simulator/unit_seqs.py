from .core import gen_unique_seq, mutate_seq


def gen_units(master_unit, edit_weights_list):
    """Generate a set of error-free unit sequences, which are `100 * p_diversity` % different
    from `master_unit`."""
    return [mutate_seq(master_unit, edit_weights) for edit_weights_list]


def 
