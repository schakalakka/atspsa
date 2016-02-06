import sys

sys.path.insert(0, '/home/andreas/workspace/swalign/')

from ext_swalign import *


def read_align_file(filename):
    return 0


def write_align_file(filename, scores):
    with open(filename + '.align', 'w') as f:
        for key, elem in scores.items():
            f.write('{}\t{}\t{}\n'.format(key[0], key[1], elem))


def compute_scores(reads, filename):
    temp_scores = {(i, j): (fast_smith_waterman(x, y)[2:5]) for i, x in enumerate(reads) for j, y in enumerate(reads)}
    with open(filename + '_alignstat.stat', 'w') as f:
        [f.write('{}\t{}\n'.format(val[1], val[2])) for key, val in temp_scores.items()]

    return {key: val[0] for key, val in temp_scores.items()}
    # scores = {}
    # for i, x in enumerate(reads):
    #     for j, y in enumerate(reads):
    #         if i < j:
    #             print(i, j)
    #             #scores[(i, j)] = align.Align(x, y, init=True).get_score()
    #             scores[(i, j)] = swalign.fast_smith_waterman(x, y)[2]
    # return scores
