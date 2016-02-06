import os
from multiprocessing import Pool

from config import *


def run_lkh(lkh_par_file):
    os.system("{}LKH {}".format(LKH_EXE, lkh_par_file))


def pool_func(current_elem):
    curr_DIR = DIR + current_elem
    os.chdir(curr_DIR)
    filename = current_elem
    run_lkh(filename + ".par")
    # run_lkh(filename + "_sparse.par")


dirs = os.listdir(DIR)

pool = Pool(processes=CPUS)
pool.map(pool_func, dirs)
