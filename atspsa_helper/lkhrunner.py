import os

from config import *


def run_lkh(lkh_par_file):
    os.system("{}LKH {} >> {}".format(LKH_EXE, lkh_par_file, lkh_par_file.split('.par')[0] + '.out'))


def pool_func(file):
    curr_DIR = DIR + file.split('/')[-2]
    os.chdir(curr_DIR)
    if not os.path.isfile(file.split('.par')[0] + '.tour') or not os.path.isfile(file.split('.par')[0] + '.out'):
        print(file)
        run_lkh(file)


# parameter_files = sorted(glob.glob(DIR+'*/*.par'))
# pool = Pool(processes=CPUS)
# pool.map(pool_func, parameter_files)
pool_func('/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign_full.par')
