import glob
import os
from multiprocessing import Pool

from config import *


def run_lkh(lkh_par_file):
    os.system("{}LKH {} > {}".format(LKH_EXE, lkh_par_file, lkh_par_file.split('.par')[0] + '.out'))


def pool_func(file):
    curr_DIR = DIR + file.split('/')[-2]
    os.chdir(curr_DIR)
    print(file)
    run_lkh(file)


if __name__ == '__main__':
    # parameter_files = sorted(glob.glob(DIR+'*/*.par'))
    # if not os.path.isfile(file.split('.par')[0] + '.tour') or not os.path.isfile(file.split('.par')[0] + '.out'):
    # parameter_files = ['/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c20_l700/calign_full.par', '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign_full.par']
    parameter_files = sorted(glob.glob(DIR + 'ref1shuffled_c[5]*/calign_full.par'))
    pool = Pool(processes=CPUS)
    pool.map(pool_func, parameter_files)
