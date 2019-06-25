from logzero import logger
from BITS.util.proc import run_command

dir_name = 'datander'


def calc_n_blocks(db_fname):
    '''Extract the number of blocks from the db file.'''
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith('blocks'):
                return int(line.split('=')[1].strip())
    logger.error(f'No information on the number of blocks in {db_fname}')


def run_datander(db_prefix, n_core=1, scheduler=None):
    run_command(f'rm -f .{db_prefix}.*.tan.* .{db_prefix}.tan.* TAN.*')
    run_command(f'mkdir -p {dir_name}; rm -f {dir_name}/*')

    n_blocks = calc_n_blocks(f'{db_prefix}.db')
    script = '\n'.join([run_command(f'HPC.TANmask -T{n_core} {db_prefix}.db'),
                        f'Catrack -v {db_prefix} tan' if n_blocks > 1 else '',
                        f'rm .{db_prefix}.*.tan.*' if n_blocks > 1 else ''])
    script_fname = f'{dir_name}/run_datander.sh'
    if scheduler is None:
        with open(script_fname, 'w') as f:
            f.write(f'{script}\n')
        run_command(f'bash {script_fname} > {dir_name}/log 2>&1')
    else:
        scheduler.submit(script,
                         script_fname,
                         job_name='datander',
                         log_fname=f'{dir_name}/log',
                         n_core=n_core,
                         wait=True)
