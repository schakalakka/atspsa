import glob

from atspsacore import createtsp, read_score_file
from config import *


def create_full_atsp(file, circular=False):
    current_file = file.split('.score')[0] + '_full'
    print(current_file)
    nr_of_reads, scores = read_score_file.read_score_file(file)
    scores[scores < MINIMAL_OVERLAP_SCORE] = -BIG_M_WEIGHT
    createtsp.write_full_atsp(current_file, nr_of_reads, scores, circular)


def create_full_atsp_via_scores(outputfile, scores, circular=False):
    nr_of_reads = len(scores)
    scores[scores < MINIMAL_OVERLAP_SCORE] = -BIG_M_WEIGHT
    createtsp.write_full_atsp(outputfile, nr_of_reads, scores, circular)


if __name__ == '__main__':
    # for file in glob.glob(DIR+'*/*.score'):
    #     create_full_atsp(file)
    for file in glob.glob(DIR + 'ref1shuffled_c[5]*/calign.score'):
        create_full_atsp(file, True)
        # create_full_atsp('/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c20_l700/calign.score', True)
        # create_full_atsp('/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.score', True)
