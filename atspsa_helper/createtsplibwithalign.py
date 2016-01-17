import os
from multiprocessing import Pool

from atspsacore import computealignment
from atspsacore import createtsp
from atspsacore import parser
from config import *


def filter_score(scores, n):
    return scores


dirs = os.listdir(DIR)


def pool_func(filename):
    curr_DIR = DIR + filename
    os.chdir(curr_DIR)
    reads = parser.parse_fasta(filename)
    scores = computealignment.compute_scores(reads)
    computealignment.write_align_file(filename, scores)
    createtsp.prepare_lkh(filename, len(reads), scores)

    # sparse_scores = {}
    # for i, elem in enumerate(reads):
    #     curr_scores = {(x, y): scores[(x, y)][0] for (x, y) in scores.keys() if i == x}
    #     vals = sorted(list(curr_scores.values()), reverse=True)[0:sparsification]
    #     sparse_scores.update({key: val for key, val in curr_scores.items() if val in vals})
    #
    # createtsp.prepare_lkh(filename + '_sparse', len(reads),
    #                       sparse_scores)


pool = Pool(processes=CPUS)
pool.map(pool_func, dirs)
