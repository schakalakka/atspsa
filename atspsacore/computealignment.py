import pickle
from ext_calign import *

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
    if ALIGNMENT_TYPE == 'semiglobal':
        temp_scores = {}
        for i, a in enumerate(reads):
            for j, b in enumerate(reads):
                if i < j:
                    (score, row, col) = align(a, b, ALIGNMENT_TYPE)[2:5]
                    alen = len(a)
                    blen = len(b)
                    if score >= MIN_SCORE:
                        if alen == row and blen == col:
                            temp_scores[(i, j)] = (score, row, col)
                            temp_scores[(j, i)] = (score, row, col)
                        elif alen == row and blen > col:
                            temp_scores[(i, j)] = (score, row, col)
                        elif alen > row and blen == col:
                            temp_scores[(j, i)] = (score, row, col)
                        else:
                            print('Huh?')
    else:
        #compute alignment scores
        #(index of read 1, index of read 2): (score, matrix row, matrix col)
        temp_scores = {(i, j): (align(x, y, ALIGNMENT_TYPE)[2:5]) for i, x in enumerate(reads) for j, y in enumerate(reads) if i < j}
        #discard all scores under threshold MIN_SCORE
        temp_scores = {key: val for key, val in temp_scores.items() if val[0] >= MIN_SCORE}
        #make scores directed according to the matrix index information




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
