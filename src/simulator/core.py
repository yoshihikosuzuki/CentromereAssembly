import random
random.seed(111)

nucleotides = ('a', 'c', 'g', 't')
edit_ops = ('=', 'X', 'I', 'D')


def set_seed(seed):
    random.seed(seed)


def gen_unique_seq(length):
    """Generate a NON tandem repeat sequence of `length` bp,
    which can be a unit sequence of a tandem repeat.
    """
    return ''.join(random.choices(nucleotides, k=length))


def gen_edit_script(length, edit_weights):
    """Generate a random edit script whose query length is `length`
    based on the weights of the edit operations `edit_weights`."""
    s = ""
    query_len = 0   # length of the query string that `s` accepts
    while query_len < length:
        c = random.choices(edit_ops, weights=edit_weights)[0]
        s += c
        if c in set(['=', 'X', 'D']):
            query_len += 1
    return s


def apply_edit_script(true_seq, edit_script):
    """Mutate `true_seq` based on `edit_script`."""
    pos = 0   # on `true_seq`
    obs_seq = ""   # mutated sequence
    for edit_op in edit_script:
        if edit_op == '=':
            obs_seq += true_seq[pos]
        elif edit_op == 'X':
            obs_seq += random.choice(list(filter(lambda n: n != true_seq[pos], nucleotides)))
        elif edit_op == 'I':
            obs_seq += random.choice(nucleotides)

        if edit_op in set(['=', 'X', 'D']):
            pos += 1

    assert pos == len(true_seq), "Given edit script does not accept the sequence"
    return obs_seq


def mutate_seq(true_seq, edit_weights):
    """Insert random mutations to `true_seq` given a list of weights (= probabilities)
    for each edit operation, i.e., {=, X, I, D} in this order."""
    return apply_edit_script(true_seq, gen_edit_script(len(true_seq), edit_weights))
