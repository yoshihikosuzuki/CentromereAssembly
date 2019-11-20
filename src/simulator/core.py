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


def gen_stochastic_edit_script(length, edit_weights):
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


def gen_deterministic_edit_script(length, edit_weights):
    """Generate an edit script which has each edit operation whose number is exactly proportional
    to `edit_weights`."""
    # Calculate the number of each edit operation to be inserted
    weight_sum = sum(edit_weights)
    edit_nums = [int(length * weight / weight_sum) for weight in edit_weights[1:]]   # for X, I, D
    # Determine the positions of the edit operations to be inserted
    s = ('X' * edit_nums[0]
         + 'I' * edit_nums[1]
         + 'D' * edit_nums[2]
         + '=' * (length - edit_nums[0] - edit_nums[2]))
    return ''.join(random.sample(s, len(s)))


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


def sequence_seq(true_seq, edit_weights):
    """Insert stochastic errors to `true_seq` given a list of weights (= probabilities)
    for each edit operation, i.e., {=, X, I, D} in this order."""
    return apply_edit_script(true_seq, gen_stochastic_edit_script(len(true_seq), edit_weights))


def mutate_seq(true_seq, edit_weights):
    """Insert deterministic mutations to `true_seq` given a list of weights (= proportions)
    into random positions."""
    return apply_edit_script(true_seq, gen_deterministic_edit_script(len(true_seq), edit_weights))
