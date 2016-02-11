import pickle
from ext_swalign import *

from config import *


def read_align_file(filename):
    return 0


def write_align_file(filename, scores):
    # writescores as pickle for easier loading
    with open(filename + '.alignpickle', 'wb') as f:
        pickle.dump(scores, f)
    # write scores as human readable file
    # tail head weight
    with open(filename + '.align', 'w') as f:
        for key, elem in scores.items():
            f.write('{}\t{}\t{}\n'.format(key[0], key[1], elem))


def compute_scores(reads, filename):
    temp_scores = {(i, j): (fast_smith_waterman(x, y)[2:5]) for i, x in enumerate(reads) for j, y in enumerate(reads)}
    temp_scores = {key: val for key, val in temp_scores.items() if key[0] > key[1] and val[0] >= MIN_SCORE}

    pre_alignstat = [(val[1], val[2]) for key, val in temp_scores.items()]
    alignstat = {elem: pre_alignstat.count(elem) for elem in set(pre_alignstat)}

    with open(filename + '_alignstat.pickle', 'wb') as f:
        pickle.dump(alignstat, f)

    return {key: val[0] for key, val in temp_scores.items()}
    # scores = {}
    # for i, x in enumerate(reads):
    #     for j, y in enumerate(reads):
    #         if i < j:
    #             print(i, j)
    #             #scores[(i, j)] = align.Align(x, y, init=True).get_score()
    #             scores[(i, j)] = swalign.fast_smith_waterman(x, y)[2]
    # return scores
